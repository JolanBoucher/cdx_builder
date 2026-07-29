//
// Created by Jolan on 2026-07-24.
//
#include "graph_linearization.h"
#include <stdexcept>
#include <string>
#include <utility>

CSRMatrix::CSRMatrix(
    std::vector<Offset> forward_key,
    std::vector<NodeId> forward_nodes,
    std::vector<Weight> forward_weights,
    std::vector<Alpha>  forward_alpha,
    std::vector<Offset> backward_key,
    std::vector<NodeId> backward_nodes,
    std::vector<Weight> backward_weights,
    std::vector<Alpha>  backward_alpha
)
    : forward_key_(std::move(forward_key)),
      forward_nodes_(std::move(forward_nodes)),
      forward_weights_(std::move(forward_weights)),
      forward_alpha_(std::move(forward_alpha)),
      backward_key_(std::move(backward_key)),
      backward_nodes_(std::move(backward_nodes)),
      backward_weights_(std::move(backward_weights)),
      backward_alpha_(std::move(backward_alpha))
{
    // Ensure all internal invariants, monotonicities, and array bounds hold upon construction
    validate();
}

// --- Fast Unchecked Edge Range ---
CSRMatrix::EdgeRange CSRMatrix::forward_edges_unchecked(const size_t nid) const noexcept {
    // Extract slice offsets directly from the prefix sum key vector
    const Offset start = forward_key_[nid];
    const Offset end   = forward_key_[nid + 1];
    const size_t count = end - start;

    if (count == 0) {
        return {nullptr, nullptr, 0};
    }

    // Return contiguous memory pointers to avoid dynamic allocations in hot relaxation loops
    return {forward_nodes_.data() + start, forward_weights_.data() + start,count
    };
}

CSRMatrix::EdgeRange CSRMatrix::backward_edges_unchecked(const size_t nid) const noexcept {
    const Offset start = backward_key_[nid];
    const Offset end   = backward_key_[nid + 1];
    const size_t count = end - start;

    if (count == 0) {
        return {nullptr, nullptr, 0};
    }

    return {backward_nodes_.data() + start, backward_weights_.data() + start,count
    };
}

// --- Checked Accessors ---
CSRMatrix::Offset CSRMatrix::forward_begin(const size_t nid) const {
    check_node_index(nid);
    return forward_key_[nid];
}

CSRMatrix::Offset CSRMatrix::forward_end(const size_t nid) const {
    check_node_index(nid);
    return forward_key_[nid + 1];
}

CSRMatrix::Offset CSRMatrix::backward_begin(const size_t nid) const {
    check_node_index(nid);
    return backward_key_[nid];
}

CSRMatrix::Offset CSRMatrix::backward_end(const size_t nid) const {
    check_node_index(nid);
    return backward_key_[nid + 1];
}

size_t CSRMatrix::forward_degree(const size_t nid) const {
    return forward_end(nid) - forward_begin(nid);
}

size_t CSRMatrix::backward_degree(const size_t nid) const {
    return backward_end(nid) - backward_begin(nid);
}

bool CSRMatrix::has_forward_edges(const size_t nid) const {
    return forward_begin(nid) != forward_end(nid);
}

bool CSRMatrix::has_backward_edges(const size_t nid) const {
    return backward_begin(nid) != backward_end(nid);
}

CSRMatrix::NodeId CSRMatrix::forward_node_at(const Offset position) const {
    check_forward_position(position);
    return forward_nodes_[static_cast<size_t>(position)];
}

CSRMatrix::Weight CSRMatrix::forward_weight_at(const Offset position) const {
    check_forward_position(position);
    return forward_weights_[static_cast<size_t>(position)];
}

CSRMatrix::NodeId CSRMatrix::backward_node_at(const Offset position) const {
    check_backward_position(position);
    return backward_nodes_[static_cast<size_t>(position)];
}

CSRMatrix::Weight CSRMatrix::backward_weight_at(const Offset position) const {
    check_backward_position(position);
    return backward_weights_[static_cast<size_t>(position)];
}

CSRMatrix::Alpha CSRMatrix::forward_alpha(const size_t nid) const {
    check_node_index(nid);
    return forward_alpha_[nid];
}

CSRMatrix::Alpha CSRMatrix::backward_alpha(const size_t nid) const {
    check_node_index(nid);
    return backward_alpha_[nid];
}

// --- Validation Logic ---
void CSRMatrix::check_node_index(const size_t nid) const {
    if (!valid_node_index(nid)) {
        throw std::out_of_range("CSRMatrix: node index " + std::to_string(nid) + " out of bounds");
    }
}

void CSRMatrix::check_forward_position(const Offset position) const {
    if (position >= forward_nodes_.size()) {
        throw std::out_of_range("CSRMatrix: forward position out of bounds");
    }
}

void CSRMatrix::check_backward_position(const Offset position) const {
    if (position >= backward_nodes_.size()) {
        throw std::out_of_range("CSRMatrix: backward position out of bounds");
    }
}

void CSRMatrix::validate() const {
    // 1. Validate key vector structure and offset integrity
    validate_key_structure();

    const size_t capacity = node_capacity();

    // 2. Validate element count alignments across forward and backward structures
    if (forward_nodes_.size() != forward_weights_.size() ||
        backward_nodes_.size() != backward_weights_.size()) {
        throw std::invalid_argument("CSRMatrix: node and weight array sizes do not match");
    }

    if (forward_nodes_.size() != backward_nodes_.size()) {
        throw std::invalid_argument("CSRMatrix: forward and backward edge counts differ");
    }

    // 3. Single optimized pass per direction to check node IDs, weights, and total sum conservation
    const uint64_t fwd_sum = validate_and_sum_forward_view(capacity);
    const uint64_t bwd_sum = validate_and_sum_backward_view(capacity);

    if (fwd_sum != bwd_sum) {
        throw std::invalid_argument("CSRMatrix: forward and backward weight sums do not match");
    }

    // 4. Ensure trust alpha factors are well-formed and bounded in [0.0, 1.0]
    validate_alphas(capacity);
}

void CSRMatrix::validate_key_structure() const {
    if (forward_key_.empty() || backward_key_.empty()) {
        throw std::invalid_argument("CSRMatrix: key arrays must not be empty");
    }

    if (forward_key_.size() != backward_key_.size()) {
        throw std::invalid_argument("CSRMatrix: forward and backward key sizes differ");
    }

    // CSR key arrays must always start with an offset of 0
    if (forward_key_.front() != 0 || backward_key_.front() != 0) {
        throw std::invalid_argument("CSRMatrix: key arrays must start at offset 0");
    }

    // Enforce strict monotonicity of CSR offset vectors
    for (size_t i = 1; i < forward_key_.size(); ++i) {
        if (forward_key_[i] < forward_key_[i - 1] || backward_key_[i] < backward_key_[i - 1]) {
            throw std::invalid_argument("CSRMatrix: key offsets are not monotonic");
        }
    }

    // The last offset in key must equal total element count in adjacent node vectors
    if (forward_key_.back() != forward_nodes_.size() || backward_key_.back() != backward_nodes_.size()) {
        throw std::invalid_argument("CSRMatrix: final key offset does not match edge array size");
    }
}

uint64_t CSRMatrix::validate_and_sum_forward_view(const size_t capacity) const {
    uint64_t sum = 0;
    for (size_t i = 0; i < forward_nodes_.size(); ++i) {
        const NodeId nid = forward_nodes_[i];
        const Weight weight = forward_weights_[i];

        // Guard against unseen sentinels or out-of-bounds target IDs
        if (nid == NODE_UNSEEN || nid >= capacity) {
            throw std::invalid_argument("CSRMatrix: invalid NodeId in forward_nodes at index " + std::to_string(i));
        }
        // Guard against zero-weight edges
        if (weight == 0) {
            throw std::invalid_argument("CSRMatrix: zero edge weight in forward_weights at index " + std::to_string(i));
        }
        sum += weight;
    }
    return sum;
}

uint64_t CSRMatrix::validate_and_sum_backward_view(const size_t capacity) const {
    uint64_t sum = 0;
    for (size_t i = 0; i < backward_nodes_.size(); ++i) {
        const NodeId nid = backward_nodes_[i];
        const Weight weight = backward_weights_[i];

        if (nid == NODE_UNSEEN || nid >= capacity) {
            throw std::invalid_argument("CSRMatrix: invalid NodeId in backward_nodes at index " + std::to_string(i));
        }
        if (weight == 0) {
            throw std::invalid_argument("CSRMatrix: zero edge weight in backward_weights at index " + std::to_string(i));
        }
        sum += weight;
    }
    return sum;
}

void CSRMatrix::validate_alphas(const size_t capacity) const {
    if (forward_alpha_.size() != capacity || backward_alpha_.size() != capacity) {
        throw std::invalid_argument("CSRMatrix: alpha array sizes must match node capacity");
    }

    // Validate probability range constraint [0.0, 1.0] for node trust factors
    for (size_t i = 0; i < capacity; ++i) {
        if (forward_alpha_[i] < 0.0f || forward_alpha_[i] > 1.0f ||
            backward_alpha_[i] < 0.0f || backward_alpha_[i] > 1.0f) {
            throw std::invalid_argument("CSRMatrix: alpha factors must be in range [0, 1]");
        }
    }
}