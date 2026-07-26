//
// Created by Jolan on 2026-07-24.
//

#ifndef CDX_BUILDER_GRAPH_LINEARIZATION_H
#define CDX_BUILDER_GRAPH_LINEARIZATION_H

#include <limits>
#include <unordered_map>
#include <vector>

class CSRMatrix {
public:
    using NodeId = uint32_t;
    using Weight = uint16_t; // we assume one edge won't ever have more than 65'535 paths going through it
    using Offset = uint32_t;
    using Alpha  = float;

    static constexpr NodeId NODE_UNSEEN = std::numeric_limits<NodeId>::max();

    // Vue contiguë légère pour la boucle rapide de relaxation
    struct EdgeRange {
        const NodeId* nodes{nullptr};
        const Weight* weights{nullptr};
        size_t count{0};
    };

    CSRMatrix() = default;

    CSRMatrix(
        std::vector<Offset> forward_key,
        std::vector<NodeId> forward_nodes,
        std::vector<Weight> forward_weights,
        std::vector<Alpha>  forward_alpha,
        std::vector<Offset> backward_key,
        std::vector<NodeId> backward_nodes,
        std::vector<Weight> backward_weights,
        std::vector<Alpha>  backward_alpha
    );

    // --- Accesseurs inline rapides ---
    [[nodiscard]] size_t node_capacity() const noexcept {
        return forward_key_.empty() ? 0 : forward_key_.size() - 1;
    }

    [[nodiscard]] size_t edge_count() const noexcept {
        return forward_nodes_.size();
    }

    [[nodiscard]] bool empty() const noexcept {
        return edge_count() == 0;
    }

    [[nodiscard]] bool valid_node_index(size_t nid) const noexcept {
        return nid < node_capacity();
    }

    // --- Accesseurs rapides non vérifiés (Unchecked / Hot Loops) ---
    [[nodiscard]] EdgeRange forward_edges_unchecked(size_t nid) const noexcept;
    [[nodiscard]] EdgeRange backward_edges_unchecked(size_t nid) const noexcept;

    // --- Accesseurs vérifiés (Debug / Interfaçage externe) ---
    [[nodiscard]] Offset forward_begin(size_t nid) const;
    [[nodiscard]] Offset forward_end(size_t nid) const;
    [[nodiscard]] Offset backward_begin(size_t nid) const;
    [[nodiscard]] Offset backward_end(size_t nid) const;

    [[nodiscard]] size_t forward_degree(size_t nid) const;
    [[nodiscard]] size_t backward_degree(size_t nid) const;
    [[nodiscard]] bool has_forward_edges(size_t nid) const;
    [[nodiscard]] bool has_backward_edges(size_t nid) const;

    [[nodiscard]] NodeId forward_node_at(Offset position) const;
    [[nodiscard]] Weight forward_weight_at(Offset position) const;
    [[nodiscard]] NodeId backward_node_at(Offset position) const;
    [[nodiscard]] Weight backward_weight_at(Offset position) const;

    [[nodiscard]] Alpha forward_alpha(size_t nid) const;
    [[nodiscard]] Alpha backward_alpha(size_t nid) const;

    // --- Getters bruts ---
    [[nodiscard]] const std::vector<Offset>& forward_key() const noexcept { return forward_key_; }
    [[nodiscard]] const std::vector<NodeId>& forward_nodes() const noexcept { return forward_nodes_; }
    [[nodiscard]] const std::vector<Weight>& forward_weights() const noexcept { return forward_weights_; }
    [[nodiscard]] const std::vector<Alpha>&  forward_alphas() const noexcept { return forward_alpha_; }

    [[nodiscard]] const std::vector<Offset>& backward_key() const noexcept { return backward_key_; }
    [[nodiscard]] const std::vector<NodeId>& backward_nodes() const noexcept { return backward_nodes_; }
    [[nodiscard]] const std::vector<Weight>& backward_weights() const noexcept { return backward_weights_; }
    [[nodiscard]] const std::vector<Alpha>&  backward_alphas() const noexcept { return backward_alpha_; }

    // Validation complète des invariants du CSR
    void validate() const;

private:
    std::vector<Offset> forward_key_;
    std::vector<NodeId> forward_nodes_;
    std::vector<Weight> forward_weights_;
    std::vector<Alpha>  forward_alpha_;

    std::vector<Offset> backward_key_;
    std::vector<NodeId> backward_nodes_;
    std::vector<Weight> backward_weights_;
    std::vector<Alpha>  backward_alpha_;

    void check_node_index(size_t nid) const;
    void check_forward_position(Offset position) const;
    void check_backward_position(Offset position) const;

    // Subroutines de validation
    void validate_key_structure() const;
    [[nodiscard]] uint64_t validate_and_sum_forward_view(size_t capacity) const;
    [[nodiscard]] uint64_t validate_and_sum_backward_view(size_t capacity) const;
    void validate_alphas(size_t capacity) const;
};

// from the edges of each node, infer local topoplogy and stock it in a csr matrix
CSRMatrix build_csr_matrix(
    const std::unordered_map<uint64_t, uint32_t>& weights
    );

/**
 * Memory-aligned compact structure representing an active node (exactly 32 bytes).
 */
struct ActiveNode {
    uint32_t nid;               // 4 octets
    uint32_t f_degree;        // 4 octets
    uint32_t b_degree;        // 4 octets
    CSRMatrix::Alpha f_alpha; // 4 octets (float)
    CSRMatrix::Alpha b_alpha; // 4 octets (float)
    CSRMatrix::Offset f_lo;   // 8 octets
    CSRMatrix::Offset b_lo;   // 8 octets
};                            // Total = 32 octets (including 4 octets of padding)


// use the csr matrix data and the node position median
// to perform a bidirectional Laplacian relaxation of the graph topology
std::pair<std::vector<uint32_t>, size_t> relax_topology(
    const CSRMatrix& matrix,
    const std::vector<uint32_t>& n_length,
    const std::vector<uint32_t>& median_pos,
    float convergence_threshold = 0.0001f,
    size_t max_iterations = 100,
    float lambda_factor = 0.5f);

#endif //CDX_BUILDER_GRAPH_LINEARIZATION_H