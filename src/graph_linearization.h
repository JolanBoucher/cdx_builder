/**
 * @file graph_linearization.h
 * @brief Graph Linearization & Coordinate Relaxation Module.
 *
 * This module provides high-performance data structures and algorithm orchestrators
 * to linearize graph nodes into continuous coordinate spaces. It features:
 *
 * 1. **Bidirectional Compressed Sparse Row (CSR) Representation (`CSRMatrix`)**:
 *    - Memory-efficient separate forward/backward adjacency storage for rapid graph traversal.
 *    - Node-level topological trust factors (alpha) and edge path weights.
 *
 * 2. **Cache-Aligned Active Node Layout (`ActiveNode`)**:
 *    - Exactly 32-byte cache-friendly structures designed for dense multithreaded accesses
 *      during relaxation passes.
 *
 * 3. **Parallel Topological Coordinate Relaxation (`relax_topology`)**:
 *    - OpenMP-accelerated iterative solver balancing local neighbor connectivity
 *      (via weighted barycenter) with fixed anchor positions (via regularization parameter lambda).
 */

#ifndef CDX_BUILDER_GRAPH_LINEARIZATION_H
#define CDX_BUILDER_GRAPH_LINEARIZATION_H

#include <limits>
#include <unordered_map>
#include <vector>
#include <cstdint>
#include <cstddef>

/**
 * @brief Bidirectional Compressed Sparse Row (CSR) matrix representation for graph topology and edge weights.
 *
 * This class optimizes memory layout and traversal performance for graph relaxation algorithms,
 * maintaining separate forward and backward adjacency lists, edge weights, and node-level alpha trust factors.
 */
class CSRMatrix {
public:
    using NodeId = uint32_t;
    using Weight = uint16_t;    // We assume one edge won't ever have more than 65'535 paths going through it
    using Offset = uint32_t;
    using Alpha  = float;

    static constexpr NodeId NODE_UNSEEN = std::numeric_limits<NodeId>::max();

    /** @brief Lightweight contiguous view structure for fast relaxation hot-loops. */
    struct EdgeRange {
        const NodeId* nodes{nullptr};   // Pointer to adjacent node IDs
        const Weight* weights{nullptr}; // Pointer to edge weights
        size_t count{0};                // Number of edges in the range
    };

    /** @brief Default constructor for an empty CSR matrix. */
    CSRMatrix() = default;

    /**
     * @brief Constructs a bidirectional CSRMatrix with preallocated vectors and validates invariants.
     *
     * @param forward_key Prefix sum offset vector for incoming/forward connections.
     * @param forward_nodes Adjacent source node IDs for forward view.
     * @param forward_weights Weights corresponding to forward edges.
     * @param forward_alpha Node-level trust alpha factors for forward orientation.
     * @param backward_key Prefix sum offset vector for outgoing/backward connections.
     * @param backward_nodes Adjacent destination node IDs for backward view.
     * @param backward_weights Weights corresponding to backward edges.
     * @param backward_alpha Node-level trust alpha factors for backward orientation.
     */
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

    /** @brief Deallocates all internal forward and backward CSR buffers and resets capacity to zero. */
    void clear() noexcept {
        // Forward storage
        forward_key_.clear();     forward_key_.shrink_to_fit();
        forward_nodes_.clear();   forward_nodes_.shrink_to_fit();
        forward_weights_.clear(); forward_weights_.shrink_to_fit();
        forward_alpha_.clear();   forward_alpha_.shrink_to_fit();

        // Backward storage
        backward_key_.clear();     backward_key_.shrink_to_fit();
        backward_nodes_.clear();   backward_nodes_.shrink_to_fit();
        backward_weights_.clear(); backward_weights_.shrink_to_fit();
        backward_alpha_.clear();   backward_alpha_.shrink_to_fit();
    }

    // --- Fast Inline Accessors ---
    /** @brief Returns the maximum number of nodes the matrix can index. */
    [[nodiscard]] size_t node_capacity() const noexcept {
        return forward_key_.empty() ? 0 : forward_key_.size() - 1;
    }
    /** @brief Returns the total number of edges stored in the matrix. */
    [[nodiscard]] size_t edge_count() const noexcept {
        return forward_nodes_.size();
    }
    /** @brief Checks whether the CSR matrix contains zero edges. */
    [[nodiscard]] bool empty() const noexcept {
        return edge_count() == 0;
    }

    /** @brief Validates if a given node ID falls within the matrix capacity bounds. */
    [[nodiscard]] bool valid_node_index(const size_t nid) const noexcept {
        return nid < node_capacity();
    }

    // --- Unchecked Accessors (Hot Loops / High Performance) ---
    /** @brief Returns a contiguous edge range view for forward edges of a node without bounds checking. */
    [[nodiscard]] EdgeRange forward_edges_unchecked(size_t nid) const noexcept;
    /** @brief Returns a contiguous edge range view for backward edges of a node without bounds checking. */
    [[nodiscard]] EdgeRange backward_edges_unchecked(size_t nid) const noexcept;

    // --- Checked Accessors (Debug / External Interfacing) ---
    /** @brief Returns the start offset index for forward edges of a node. */
    [[nodiscard]] Offset forward_begin(size_t nid) const;
    /** @brief Returns the end offset index for forward edges of a node. */
    [[nodiscard]] Offset forward_end(size_t nid) const;
    /** @brief Returns the start offset index for backward edges of a node. */
    [[nodiscard]] Offset backward_begin(size_t nid) const;
    /** @brief Returns the end offset index for backward edges of a node. */
    [[nodiscard]] Offset backward_end(size_t nid) const;
    /** @brief Returns the forward degree (number of incoming edges) for a node. */
    [[nodiscard]] size_t forward_degree(size_t nid) const;
    /** @brief Returns the backward degree (number of outgoing edges) for a node. */
    [[nodiscard]] size_t backward_degree(size_t nid) const;
    /** @brief Checks if a node has any forward edges. */
    [[nodiscard]] bool has_forward_edges(size_t nid) const;
    /** @brief Checks if a node has any backward edges. */
    [[nodiscard]] bool has_backward_edges(size_t nid) const;
    /** @brief Retrieves the adjacent node ID at a specific forward global table position. */
    [[nodiscard]] NodeId forward_node_at(Offset position) const;
    /** @brief Retrieves the edge weight at a specific forward global table position. */
    [[nodiscard]] Weight forward_weight_at(Offset position) const;
    /** @brief Retrieves the adjacent node ID at a specific backward global table position. */
    [[nodiscard]] NodeId backward_node_at(Offset position) const;
    /** @brief Retrieves the edge weight at a specific backward global table position. */
    [[nodiscard]] Weight backward_weight_at(Offset position) const;
    /** @brief Returns the forward alpha trust factor for a specific node. */
    [[nodiscard]] Alpha forward_alpha(size_t nid) const;
    /** @brief Returns the backward alpha trust factor for a specific node. */
    [[nodiscard]] Alpha backward_alpha(size_t nid) const;

    // --- Raw Vector Getters ---
    [[nodiscard]] const std::vector<Offset>& forward_key() const noexcept { return forward_key_; }
    [[nodiscard]] const std::vector<NodeId>& forward_nodes() const noexcept { return forward_nodes_; }
    [[nodiscard]] const std::vector<Weight>& forward_weights() const noexcept { return forward_weights_; }
    [[nodiscard]] const std::vector<Alpha>&  forward_alphas() const noexcept { return forward_alpha_; }

    [[nodiscard]] const std::vector<Offset>& backward_key() const noexcept { return backward_key_; }
    [[nodiscard]] const std::vector<NodeId>& backward_nodes() const noexcept { return backward_nodes_; }
    [[nodiscard]] const std::vector<Weight>& backward_weights() const noexcept { return backward_weights_; }
    [[nodiscard]] const std::vector<Alpha>&  backward_alphas() const noexcept { return backward_alpha_; }

    /**
     * @brief Performs complete structural and boundary validation of all CSR invariants.
     * @throws std::runtime_error If any structural invariant or boundary constraint is violated.
     */
    void validate() const;

private:
    // Forward adjacency storage structures
    std::vector<Offset> forward_key_;
    std::vector<NodeId> forward_nodes_;
    std::vector<Weight> forward_weights_;
    std::vector<Alpha>  forward_alpha_;

    // Backward adjacency storage structures
    std::vector<Offset> backward_key_;
    std::vector<NodeId> backward_nodes_;
    std::vector<Weight> backward_weights_;
    std::vector<Alpha>  backward_alpha_;

    // Internal safety check routines
    void check_node_index(size_t nid) const;
    void check_forward_position(Offset position) const;
    void check_backward_position(Offset position) const;

    // Internal validation subroutines
    void validate_key_structure() const;
    [[nodiscard]] uint64_t validate_and_sum_forward_view(size_t capacity) const;
    [[nodiscard]] uint64_t validate_and_sum_backward_view(size_t capacity) const;
    void validate_alphas(size_t capacity) const;
};

/**
 * @brief Public orchestrator function for constructing the complete bidirectional linearization CSRMatrix.
 *
 * Coordinates the full multi-pass construction pipeline by:
 * 1. Building the prefix sum key offsets for forward and backward directions via `build_csr_key`.
 * 2. Populating the underlying adjacency nodes and weights buffers via `fill_csr_matrix`.
 * 3. Calculating the topological trust alpha factors for every node via `compute_alpha_factors`.
 * 4. Wrapping and validating the finalized components into the target `CSRMatrix` structure.
 *
 * @param weights Reference to the unordered map storing packed 64-bit node pairs and their corresponding weights.
 * @return CSRMatrix The fully constructed and validated bidirectional compressed sparse row matrix object.
 */
CSRMatrix build_csr_matrix(
    const std::unordered_map<uint64_t, uint32_t>& weights);

/**
 * @brief Memory-aligned compact structure representing an active node in the graph (exactly 32 bytes).
 *
 * Caches essential node metadata and memory offsets to optimize traversal and
 * coordinate updates during the topological relaxation hot-loop.
 */
struct ActiveNode {
    uint32_t nid;               // Node identifier (4 bytes)
    uint32_t f_degree;          // Forward incoming edge degree (4 bytes)
    uint32_t b_degree;          // Backward outgoing edge degree (4 bytes)
    CSRMatrix::Alpha f_alpha;   // Forward topological alpha trust factor (4 bytes, float)
    CSRMatrix::Alpha b_alpha;   // Backward topological alpha trust factor (4 bytes, float)
    CSRMatrix::Offset f_lo;     // Low offset pointer for forward adjacency slice (8 bytes)
    CSRMatrix::Offset b_lo;     // Low offset pointer for backward adjacency slice (8 bytes)
};                              // Total size = 32 bytes (includes 4 bytes of struct alignment padding)


/**
 * @brief Shifts relaxed coordinates by a uniform constant so the minimum active-node position is >= 0.
 *
 * `relax_topology`'s update rule has no explicit lower bound: each node's new position is a
 * (lambda-regularized) weighted barycenter of its neighbors, offset by +/- their sequence
 * lengths. When `lambda_factor` is small, the anchor term barely pulls positions back toward
 * their original non-negative median offsets, so the system behaves like pure Laplacian
 * smoothing -- which is invariant to a uniform translation of *all* coordinates (shifting every
 * node by the same constant doesn't change any of the barycenter differences that drive
 * convergence). The forward update adds a neighbor's length while the backward update subtracts
 * the current node's length, so this "gauge" isn't perfectly symmetric; over enough iterations
 * (as happens with a very low lambda, which converges slowly) that asymmetry lets the whole
 * system drift in one direction until part of it lands below zero.
 *
 * Since only *relative* positions carry meaning here, the gauge is fixed post-hoc instead of
 * being constrained iteration-by-iteration (which would distort local structure): this function
 * finds the minimum position among active nodes and, if it's negative, shifts every active node
 * by the same amount so the whole system sits at or above zero again. This preserves every
 * pairwise distance and ordering computed during relaxation. Nodes not listed in `active_nodes`
 * (isolated or non-existent nodes) are left untouched.
 *
 * @param active_nodes The active node metadata produced by relax_topology's node-filtering pass.
 * @param nodes_pos In/out vector of continuous node positions, indexed by node id. Only entries
 *        referenced by `active_nodes` are read or modified.
 * @return double The magnitude of the shift that was applied (0.0 if none was necessary, e.g.
 *         `active_nodes` is empty or the minimum position was already non-negative).
 */
double reanchor_to_nonnegative_origin(
    const std::vector<ActiveNode>& active_nodes,
    std::vector<double>& nodes_pos);

/**
 * @brief Executes the topological relaxation algorithm to optimize continuous node coordinates.
 *
 * Performs an iterative, parallelized coordinate relaxation process over the graph nodes
 * using a bidirectional compressed sparse row (CSR) matrix. The algorithm balances path-derived
 * anchor positions with local topological neighbor connectivity, regularized by alpha confidence
 * factors and a global lambda anchor parameter.
 *
 * Mathematically, for each active node v, forward and backward predictions are computed
 * via weighted barycenters of neighboring positions, combined with topological alphas (alpha),
 * and regularized against the fixed anchor position:
 * - x_local = alpha * barycenter + (1.0 - alpha) * x_current
 * - x_pred = (1.0 - lambda) * x_local + lambda * x_anchor
 *
 * The system tracks convergence using Root Mean Square (RMS) displacement across active nodes,
 * terminating early if stability falls below the given `convergence_threshold` or when
 * `max_iterations` is reached. OpenMP is utilized for parallel execution on large node sets.
 *
 * Once the loop finishes, the resulting coordinates are passed through
 * `reanchor_to_nonnegative_origin` (see its docstring) to correct for coordinate drift, which
 * becomes especially pronounced at low `lambda_factor` values.
 *
 * @param matrix The bidirectional `CSRMatrix` containing graph adjacencies, weights, and alpha factors.
 * @param n_length Vector mapping node indices to their lengths (unseen nodes marked with `NODE_UNSEEN_32`).
 * @param median_pos Vector of initial median coordinate positions used as anchors.
 * @param convergence_threshold Early-stopping threshold for mean RMS displacement.
 * @param max_iterations Maximum allowed relaxation iterations.
 * @param lambda_factor Regularization coefficient controlling anchor adherence in the range [0.0, 1.0].
 * @return std::pair<std::vector<uint32_t>, size_t> A pair containing the finalized 32-bit integer
 *         coordinates vector and the total number of executed iterations.
 */
std::pair<std::vector<uint32_t>, size_t> relax_topology(
    const CSRMatrix& matrix,
    const std::vector<uint32_t>& n_length,
    const std::vector<uint32_t>& median_pos,
    float convergence_threshold = 0.0001f,
    size_t max_iterations = 100,
    float lambda_factor = 0.5f);

#endif //CDX_BUILDER_GRAPH_LINEARIZATION_H
