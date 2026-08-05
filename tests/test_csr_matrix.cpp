/**
 * @file test_csr_matrix.cpp
 * @brief Unit test suite for the CSRMatrix class (src/graph_linearization.h / CSRMatrix.cpp).
 *
 * Covers, in order:
 *   - Default construction and empty-state accessors.
 *   - Happy-path bidirectional construction (valid invariants).
 *   - Every distinct invariant enforced by validate() / validate_key_structure() /
 *     validate_and_sum_forward_view() / validate_and_sum_backward_view() / validate_alphas(),
 *     one test per throw site.
 *   - Unchecked hot-loop accessors (forward_edges_unchecked / backward_edges_unchecked).
 *   - Checked accessors (begin/end/degree/has_edges/node_at/weight_at/alpha) and their
 *     out_of_range behavior on invalid input.
 *   - clear().
 *
 * This module has zero external dependencies (no gbwt/gbwtgraph/libbdsg), which is why
 * tests/CMakeLists.txt links test_csr_matrix only against CSRMatrix.cpp; build_csr_matrix()
 * and relax_topology() (also declared in graph_linearization.h, defined in
 * graph_linearization.cpp) are intentionally out of scope here and covered by
 * test_graph_linearization instead.
 *
 * Shared test infrastructure:
 *   - CsrFixtureBuilder: builds small, explicit, valid bidirectional CSRMatrix instances
 *     (as plain vectors) for a tiny 3-node graph, so individual tests can start from a known
 *     -good baseline and mutate exactly one field to trigger a specific invariant violation.
 */

#include <gtest/gtest.h>
#include <limits>
#include <vector>
#include "graph_linearization.h"

namespace {

using NodeId = CSRMatrix::NodeId;
using Weight = CSRMatrix::Weight;
using Offset = CSRMatrix::Offset;
using Alpha  = CSRMatrix::Alpha;

/**
 * @brief Bundles the eight raw CSR vectors so individual tests can copy a known-valid
 *        baseline and mutate exactly one field before constructing a CSRMatrix.
 *
 * Baseline topology (3 nodes: 0, 1, 2):
 *   Forward view (per-node incoming edges): node 1 <- node 0 (weight 5), node 2 <- node 0 (weight 3)
 *   Backward view (per-node outgoing edges): node 0 -> node 1 (weight 5), node 0 -> node 2 (weight 3)
 * This mirrors the real production convention: forward/backward describe the same edges from
 * opposite endpoints, so their per-direction weight sums must agree (5 + 3 == 5 + 3).
 */
struct CsrFixtureBuilder {
    std::vector<Offset> forward_key    = {0, 0, 1, 2};
    std::vector<NodeId> forward_nodes  = {0, 0};
    std::vector<Weight> forward_weights = {5, 3};
    std::vector<Alpha>  forward_alpha  = {1.0f, 0.5f, 0.5f};

    std::vector<Offset> backward_key    = {0, 2, 2, 2};
    std::vector<NodeId> backward_nodes  = {1, 2};
    std::vector<Weight> backward_weights = {5, 3};
    std::vector<Alpha>  backward_alpha  = {1.0f, 0.0f, 0.0f};

    [[nodiscard]] CSRMatrix build() const {
        return CSRMatrix(
            forward_key, forward_nodes, forward_weights, forward_alpha,
            backward_key, backward_nodes, backward_weights, backward_alpha
        );
    }
};

} // namespace

// ============================================================================
// 1. DEFAULT CONSTRUCTION / EMPTY STATE
// ============================================================================

TEST(CSRMatrixDefaultTest, DefaultConstructedMatrixIsEmpty) {
    const CSRMatrix matrix;
    EXPECT_EQ(matrix.node_capacity(), 0u);
    EXPECT_EQ(matrix.edge_count(), 0u);
    EXPECT_TRUE(matrix.empty());
    EXPECT_FALSE(matrix.valid_node_index(0));
}

// ============================================================================
// 2. HAPPY-PATH CONSTRUCTION
// ============================================================================

TEST(CSRMatrixConstructionTest, ValidBidirectionalMatrixConstructsAndExposesRawVectors) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    EXPECT_EQ(matrix.node_capacity(), 3u);   // forward_key.size() - 1
    EXPECT_EQ(matrix.edge_count(), 2u);      // forward_nodes.size()
    EXPECT_FALSE(matrix.empty());

    EXPECT_EQ(matrix.forward_key(), fixture.forward_key);
    EXPECT_EQ(matrix.forward_nodes(), fixture.forward_nodes);
    EXPECT_EQ(matrix.forward_weights(), fixture.forward_weights);
    EXPECT_EQ(matrix.forward_alphas(), fixture.forward_alpha);

    EXPECT_EQ(matrix.backward_key(), fixture.backward_key);
    EXPECT_EQ(matrix.backward_nodes(), fixture.backward_nodes);
    EXPECT_EQ(matrix.backward_weights(), fixture.backward_weights);
    EXPECT_EQ(matrix.backward_alphas(), fixture.backward_alpha);
}

TEST(CSRMatrixConstructionTest, MatrixWithNoEdgesButNonEmptyKeysIsValid) {
    // 2 nodes, zero edges: keys are all-zero prefix sums, node/weight arrays empty.
    const CSRMatrix matrix(
        {0, 0, 0}, {}, {}, {0.5f, 0.5f},
        {0, 0, 0}, {}, {}, {0.5f, 0.5f}
    );
    EXPECT_EQ(matrix.node_capacity(), 2u);
    EXPECT_EQ(matrix.edge_count(), 0u);
    EXPECT_TRUE(matrix.empty());
}

// ============================================================================
// 3. validate() -- ONE TEST PER ENFORCED INVARIANT
// ============================================================================
//
// Each test starts from CsrFixtureBuilder's known-valid baseline and mutates exactly
// one field so the resulting construction throws std::invalid_argument (or
// std::out_of_range where the failure is a bounds issue), matching the throw site
// documented above each test.

TEST(CSRMatrixValidateTest, EmptyForwardKeyThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_key.clear();
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, EmptyBackwardKeyThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_key.clear();
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, MismatchedKeySizesThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_key.push_back(2); // now size 5 vs forward's size 4
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardKeyNotStartingAtZeroThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_key.front() = 1;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, BackwardKeyNotStartingAtZeroThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_key.front() = 1;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, NonMonotonicForwardKeyThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_key = {0, 2, 1, 2}; // decreases at index 2
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, NonMonotonicBackwardKeyThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_key = {0, 2, 1, 2}; // decreases at index 2
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, FinalForwardKeyOffsetMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_key.back() = 3; // does not match forward_nodes.size() == 2
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, FinalBackwardKeyOffsetMismatchThrows) {
    CsrFixtureBuilder fixture;
    // Increase (not decrease) the final offset so monotonicity still holds and this
    // isolates the final-offset-vs-edge-count check rather than the monotonic one.
    fixture.backward_key.back() = 3; // does not match backward_nodes.size() == 2
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardNodesWeightsSizeMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_weights.push_back(1); // 3 weights vs 2 nodes
    // Keep the key's final offset consistent with forward_nodes so this test isolates
    // the nodes/weights mismatch rather than the key-structure check above.
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, BackwardNodesWeightsSizeMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_weights.push_back(1); // 3 weights vs 2 nodes
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardBackwardEdgeCountMismatchThrows) {
    CsrFixtureBuilder fixture;
    // Add a third forward edge (and matching weight/key) but leave backward at 2 edges.
    fixture.forward_nodes.push_back(2);
    fixture.forward_weights.push_back(1);
    fixture.forward_key = {0, 0, 1, 3};
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, NodeUnseenSentinelInForwardNodesThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_nodes[0] = CSRMatrix::NODE_UNSEEN;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, OutOfBoundsNodeIdInForwardNodesThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_nodes[0] = 99; // >= node_capacity() == 3
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, NodeUnseenSentinelInBackwardNodesThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_nodes[0] = CSRMatrix::NODE_UNSEEN;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, OutOfBoundsNodeIdInBackwardNodesThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_nodes[0] = 99; // >= node_capacity() == 3
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ZeroWeightInForwardWeightsThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_weights[0] = 0;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ZeroWeightInBackwardWeightsThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_weights[0] = 0;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardBackwardWeightSumMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_weights[0] = 6; // forward sum becomes 9, backward sum stays 8
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardAlphaSizeMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_alpha.pop_back(); // size 2 vs node_capacity() == 3
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, BackwardAlphaSizeMismatchThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_alpha.pop_back(); // size 2 vs node_capacity() == 3
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardAlphaBelowZeroThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_alpha[0] = -0.01f;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, ForwardAlphaAboveOneThrows) {
    CsrFixtureBuilder fixture;
    fixture.forward_alpha[0] = 1.01f;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, BackwardAlphaOutOfRangeThrows) {
    CsrFixtureBuilder fixture;
    fixture.backward_alpha[0] = 2.0f;
    EXPECT_THROW(fixture.build(), std::invalid_argument);
}

TEST(CSRMatrixValidateTest, BoundaryAlphaValuesZeroAndOneAreAccepted) {
    CsrFixtureBuilder fixture;
    fixture.forward_alpha = {0.0f, 1.0f, 0.0f};
    fixture.backward_alpha = {1.0f, 0.0f, 1.0f};
    EXPECT_NO_THROW(fixture.build());
}

// ============================================================================
// 4. UNCHECKED HOT-LOOP ACCESSORS (forward_edges_unchecked / backward_edges_unchecked)
// ============================================================================

TEST(CSRMatrixEdgeRangeTest, NodeWithNoForwardEdgesReturnsEmptyRange) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    // Node 0 has no *incoming* (forward) edges in the baseline fixture.
    const CSRMatrix::EdgeRange range = matrix.forward_edges_unchecked(0);
    EXPECT_EQ(range.count, 0u);
    EXPECT_EQ(range.nodes, nullptr);
    EXPECT_EQ(range.weights, nullptr);
}

TEST(CSRMatrixEdgeRangeTest, NodeWithForwardEdgesReturnsCorrectSlice) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    // Node 1 has exactly one incoming edge from node 0, weight 5 (offset 0..1).
    const CSRMatrix::EdgeRange range = matrix.forward_edges_unchecked(1);
    ASSERT_EQ(range.count, 1u);
    EXPECT_EQ(range.nodes[0], 0u);
    EXPECT_EQ(range.weights[0], 5u);
}

TEST(CSRMatrixEdgeRangeTest, NodeWithMultipleForwardEdgesReturnsFullSlice) {
    // 3 nodes; two physical edges 0->2 (weight 4) and 1->2 (weight 6), so node 2
    // has two incoming (forward) edges, matched by two outgoing (backward) edges
    // split across nodes 0 and 1 -- total edge count must agree between views.
    const CSRMatrix matrix(
        {0, 0, 0, 2}, {0, 1}, {4, 6}, {1.0f, 1.0f, 1.0f},
        {0, 1, 2, 2}, {2, 2}, {4, 6}, {1.0f, 1.0f, 1.0f}
    );
    const CSRMatrix::EdgeRange range = matrix.forward_edges_unchecked(2);
    ASSERT_EQ(range.count, 2u);
    EXPECT_EQ(range.nodes[0], 0u);
    EXPECT_EQ(range.weights[0], 4u);
    EXPECT_EQ(range.nodes[1], 1u);
    EXPECT_EQ(range.weights[1], 6u);
}

TEST(CSRMatrixEdgeRangeTest, NodeWithNoBackwardEdgesReturnsEmptyRange) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    // Node 1 has no *outgoing* (backward) edges in the baseline fixture.
    const CSRMatrix::EdgeRange range = matrix.backward_edges_unchecked(1);
    EXPECT_EQ(range.count, 0u);
    EXPECT_EQ(range.nodes, nullptr);
    EXPECT_EQ(range.weights, nullptr);
}

TEST(CSRMatrixEdgeRangeTest, NodeWithBackwardEdgesReturnsCorrectSlice) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    // Node 0 has two outgoing edges: to node 1 (weight 5) and node 2 (weight 3).
    const CSRMatrix::EdgeRange range = matrix.backward_edges_unchecked(0);
    ASSERT_EQ(range.count, 2u);
    EXPECT_EQ(range.nodes[0], 1u);
    EXPECT_EQ(range.weights[0], 5u);
    EXPECT_EQ(range.nodes[1], 2u);
    EXPECT_EQ(range.weights[1], 3u);
}

// ============================================================================
// 5. CHECKED ACCESSORS (begin/end/degree/has_edges/node_at/weight_at/alpha)
// ============================================================================

TEST(CSRMatrixCheckedAccessorTest, BeginEndAndDegreeMatchBaselineTopology) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    EXPECT_EQ(matrix.forward_begin(1), 0u);
    EXPECT_EQ(matrix.forward_end(1), 1u);
    EXPECT_EQ(matrix.forward_degree(1), 1u);
    EXPECT_TRUE(matrix.has_forward_edges(1));

    EXPECT_EQ(matrix.forward_begin(0), 0u);
    EXPECT_EQ(matrix.forward_end(0), 0u);
    EXPECT_EQ(matrix.forward_degree(0), 0u);
    EXPECT_FALSE(matrix.has_forward_edges(0));

    EXPECT_EQ(matrix.backward_begin(0), 0u);
    EXPECT_EQ(matrix.backward_end(0), 2u);
    EXPECT_EQ(matrix.backward_degree(0), 2u);
    EXPECT_TRUE(matrix.has_backward_edges(0));

    EXPECT_EQ(matrix.backward_begin(1), 2u);
    EXPECT_EQ(matrix.backward_end(1), 2u);
    EXPECT_EQ(matrix.backward_degree(1), 0u);
    EXPECT_FALSE(matrix.has_backward_edges(1));
}

TEST(CSRMatrixCheckedAccessorTest, NodeAndWeightAtReturnCorrectValues) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    EXPECT_EQ(matrix.forward_node_at(0), 0u);
    EXPECT_EQ(matrix.forward_weight_at(0), 5u);
    EXPECT_EQ(matrix.forward_node_at(1), 0u);
    EXPECT_EQ(matrix.forward_weight_at(1), 3u);

    EXPECT_EQ(matrix.backward_node_at(0), 1u);
    EXPECT_EQ(matrix.backward_weight_at(0), 5u);
    EXPECT_EQ(matrix.backward_node_at(1), 2u);
    EXPECT_EQ(matrix.backward_weight_at(1), 3u);
}

TEST(CSRMatrixCheckedAccessorTest, ForwardAndBackwardAlphaReturnCorrectValues) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    EXPECT_FLOAT_EQ(matrix.forward_alpha(0), 1.0f);
    EXPECT_FLOAT_EQ(matrix.forward_alpha(1), 0.5f);
    EXPECT_FLOAT_EQ(matrix.backward_alpha(0), 1.0f);
    EXPECT_FLOAT_EQ(matrix.backward_alpha(2), 0.0f);
}

TEST(CSRMatrixCheckedAccessorTest, ValidNodeIndexBoundaryIsExact) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();

    EXPECT_TRUE(matrix.valid_node_index(2));  // last valid index (capacity 3)
    EXPECT_FALSE(matrix.valid_node_index(3)); // one past the end
}

TEST(CSRMatrixCheckedAccessorTest, OutOfRangeNodeIndexThrowsOnEveryCheckedAccessor) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();
    constexpr size_t out_of_range_nid = 3; // node_capacity() == 3, so valid indices are [0,2]

    EXPECT_THROW(matrix.forward_begin(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.forward_end(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.backward_begin(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.backward_end(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.forward_degree(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.backward_degree(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.has_forward_edges(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.has_backward_edges(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.forward_alpha(out_of_range_nid), std::out_of_range);
    EXPECT_THROW(matrix.backward_alpha(out_of_range_nid), std::out_of_range);
}

TEST(CSRMatrixCheckedAccessorTest, OutOfRangePositionThrowsOnNodeAndWeightAt) {
    const CsrFixtureBuilder fixture;
    const CSRMatrix matrix = fixture.build();
    constexpr Offset out_of_range_position = 2; // forward/backward_nodes.size() == 2

    EXPECT_THROW(matrix.forward_node_at(out_of_range_position), std::out_of_range);
    EXPECT_THROW(matrix.forward_weight_at(out_of_range_position), std::out_of_range);
    EXPECT_THROW(matrix.backward_node_at(out_of_range_position), std::out_of_range);
    EXPECT_THROW(matrix.backward_weight_at(out_of_range_position), std::out_of_range);
}

// ============================================================================
// 6. clear()
// ============================================================================

TEST(CSRMatrixClearTest, ClearResetsCapacityAndEdgeCountToZero) {
    const CsrFixtureBuilder fixture;
    CSRMatrix matrix = fixture.build();
    ASSERT_GT(matrix.node_capacity(), 0u);
    ASSERT_GT(matrix.edge_count(), 0u);

    matrix.clear();

    EXPECT_EQ(matrix.node_capacity(), 0u);
    EXPECT_EQ(matrix.edge_count(), 0u);
    EXPECT_TRUE(matrix.empty());
    EXPECT_TRUE(matrix.forward_key().empty());
    EXPECT_TRUE(matrix.forward_nodes().empty());
    EXPECT_TRUE(matrix.forward_weights().empty());
    EXPECT_TRUE(matrix.forward_alphas().empty());
    EXPECT_TRUE(matrix.backward_key().empty());
    EXPECT_TRUE(matrix.backward_nodes().empty());
    EXPECT_TRUE(matrix.backward_weights().empty());
    EXPECT_TRUE(matrix.backward_alphas().empty());
}
