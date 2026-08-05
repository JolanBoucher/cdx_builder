/**
 * @file test_graph_linearization.cpp
 * @brief Unit test suite for the graph linearization pipeline (src/graph_linearization.h/.cpp).
 *
 * Covers, in order:
 *   - build_csr_matrix (CSR construction + alpha trust-factor computation from an edge-weight map)
 *   - relax_topology   (iterative bidirectional coordinate relaxation solver)
 *
 * `build_csr_key`, `fill_csr_matrix`, `compute_alpha_factors` and `double2uint32` are file-local
 * (declared and defined only in graph_linearization.cpp, not exposed in the header), so they are
 * exercised indirectly here through `build_csr_matrix()` / `relax_topology()`'s public contract
 * rather than called directly. CSRMatrix's own invariants (constructor validation, accessors) are
 * covered by test_csr_matrix.cpp and are assumed correct here.
 *
 * UNDOCUMENTED-BY-THE-TYPE-SYSTEM PRECONDITION (not exercised here): `build_csr_key()` guards
 * against node ids >= cfg::ARRAY_SIZE (`if (v < cfg::ARRAY_SIZE) forward_key[v]++;`), but
 * `fill_csr_matrix()` has no equivalent guard -- it unconditionally does
 * `forward_ptr[v]++`/`backward_ptr[u]++` for every entry in `weights`. In isolation this would be
 * an out-of-bounds write for an out-of-range node id. In practice this is safe because
 * cfg::ARRAY_SIZE is always set to (max node id in the loaded GBZ + 1) before any of this pipeline
 * runs, so every node id appearing in a `weights` map derived from that same graph is guaranteed to
 * be < cfg::ARRAY_SIZE by construction. Not tested with an out-of-range key here since doing so
 * would violate that caller-side invariant rather than exercise a real code path.
 *
 * Shared test infrastructure:
 *   - CfgFixture: saves/restores cfg::ARRAY_SIZE and cfg::N_HAPLO around each test (both are read
 *     by build_csr_matrix()/relax_topology(); see test_gbz_io.cpp for the rationale).
 */

#include <gtest/gtest.h>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <vector>
#include "graph_linearization.h"
#include "constant.h"

namespace {

/**
 * @brief Packs a (source, destination) 32-bit node id pair into the 64-bit edge key format
 *        expected by build_csr_matrix()'s `weights` parameter (src in the high 32 bits, dst in
 *        the low 32 bits). This is a local re-implementation of gbz_IO.h's
 *        `detail::pack_node_pair` -- deliberately not including gbz_IO.h here, since that would
 *        pull in the gbwt/gbwtgraph dependency this test target is designed to avoid (see
 *        tests/CMakeLists.txt: test_graph_linearization links neither gbwt_lib nor gbwtgraph_lib).
 */
constexpr uint64_t pack_edge(const uint32_t src, const uint32_t dst) {
    return (static_cast<uint64_t>(src) << 32) | static_cast<uint64_t>(dst);
}

/// Saves/restores the cfg:: globals this module reads (cfg::ARRAY_SIZE, cfg::N_HAPLO).
class CfgFixture : public ::testing::Test {
protected:
    void SetUp() override {
        saved_array_size_ = cfg::ARRAY_SIZE;
        saved_n_haplo_ = cfg::N_HAPLO;
    }

    void TearDown() override {
        cfg::ARRAY_SIZE = saved_array_size_;
        cfg::N_HAPLO = saved_n_haplo_;
    }

private:
    size_t saved_array_size_{};
    size_t saved_n_haplo_{};
};

} // namespace

// ============================================================================
// 1. CONSTRUCTION DU CSR + FACTEURS ALPHA (build_csr_matrix)
// ============================================================================
//
// Invariants / contract:
//   - Consumes a packed-edge-key -> weight map (see the local pack_edge() helper above) and produces a fully
//     validated bidirectional CSRMatrix (construction runs CSRMatrix::validate() internally, so a
//     malformed result would already throw before any test assertion runs).
//   - Per-node alpha trust factor formula (see graph_linearization.cpp doc comment):
//       D (dominance)  = max_edge_weight / sum_edge_weights
//       C (complexity) = max(0, (degree - N_HAPLO) / degree)
//       alpha = clamp(0.50 + 0.45 * D * (1 - C), 0.50, 0.95)
//     with alpha defaulting to 0.95 for any node with zero degree in that direction (source/sink),
//     and node index 0 always keeping the default 0.95 in both directions since the computation
//     loop starts at v = 1 (index 0 is the reserved/unused node id in this codebase's convention).

TEST_F(CfgFixture, EmptyWeightsProducesEmptyMatrixOfDeclaredCapacity) {
    cfg::ARRAY_SIZE = 3; // nodes exist (e.g. isolated), but no edges at all
    cfg::N_HAPLO = 1;

    const std::unordered_map<uint64_t, uint32_t> weights;
    const CSRMatrix matrix = build_csr_matrix(weights);

    EXPECT_EQ(matrix.node_capacity(), 3u);
    EXPECT_EQ(matrix.edge_count(), 0u);
    EXPECT_TRUE(matrix.empty());
}

TEST_F(CfgFixture, SingleEdgeProducesMatchingForwardAndBackwardViews) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 7;

    const CSRMatrix matrix = build_csr_matrix(weights);

    ASSERT_EQ(matrix.edge_count(), 1u);
    // Forward view is indexed by destination: node 2 has one incoming edge from node 1.
    ASSERT_EQ(matrix.forward_degree(2), 1u);
    EXPECT_EQ(matrix.forward_node_at(matrix.forward_begin(2)), 1u);
    EXPECT_EQ(matrix.forward_weight_at(matrix.forward_begin(2)), 7u);
    // Backward view is indexed by source: node 1 has one outgoing edge to node 2.
    ASSERT_EQ(matrix.backward_degree(1), 1u);
    EXPECT_EQ(matrix.backward_node_at(matrix.backward_begin(1)), 2u);
    EXPECT_EQ(matrix.backward_weight_at(matrix.backward_begin(1)), 7u);
}

TEST_F(CfgFixture, FanInAndFanOutProduceCorrectDegrees) {
    cfg::ARRAY_SIZE = 4;
    cfg::N_HAPLO = 2;

    // Bubble: 1 -> 3, 2 -> 3 (fan-in on node 3); node 1 also -> nothing else (no fan-out here).
    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 3)] = 2;
    weights[pack_edge(2, 3)] = 5;

    const CSRMatrix matrix = build_csr_matrix(weights);

    EXPECT_EQ(matrix.forward_degree(3), 2u); // fan-in: node 3 sees both node 1 and node 2
    EXPECT_EQ(matrix.backward_degree(1), 1u);
    EXPECT_EQ(matrix.backward_degree(2), 1u);
}

TEST_F(CfgFixture, DominantSingleEdgeYieldsMaximalAlpha) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 5; // degree (1) well below N_HAPLO, so complexity penalty C == 0

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 10; // only incoming edge to node 2: dominance D == 1.0

    const CSRMatrix matrix = build_csr_matrix(weights);

    // alpha = clamp(0.50 + 0.45 * 1.0 * (1 - 0), 0.50, 0.95) == 0.95
    EXPECT_FLOAT_EQ(matrix.forward_alpha(2), 0.95f);
}

TEST_F(CfgFixture, HighDegreeBeyondHaplotypeCountPenalizesAlpha) {
    cfg::ARRAY_SIZE = 5;
    cfg::N_HAPLO = 1; // degree (3) exceeds N_HAPLO, triggering the complexity penalty

    // Three equal-weight incoming edges to node 4 (from nodes 1, 2, 3): D = 1/3, C = (3-1)/3 = 2/3.
    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 4)] = 1;
    weights[pack_edge(2, 4)] = 1;
    weights[pack_edge(3, 4)] = 1;

    const CSRMatrix matrix = build_csr_matrix(weights);

    // alpha = clamp(0.50 + 0.45 * (1/3) * (1 - 2/3), 0.50, 0.95)
    //       = clamp(0.50 + 0.45 * (1/3) * (1/3), ...) = 0.50 + 0.05 = 0.55
    EXPECT_NEAR(matrix.forward_alpha(4), 0.55f, 1e-4f);
}

TEST_F(CfgFixture, ZeroDegreeNodeDefaultsToMaximalAlphaInThatDirection) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 4;

    const CSRMatrix matrix = build_csr_matrix(weights);

    // Node 1 has no *incoming* edges (source node): forward alpha defaults to 0.95.
    EXPECT_FLOAT_EQ(matrix.forward_alpha(1), 0.95f);
    // Node 2 has no *outgoing* edges (sink node): backward alpha defaults to 0.95.
    EXPECT_FLOAT_EQ(matrix.backward_alpha(2), 0.95f);
}

TEST_F(CfgFixture, NodeZeroAlwaysKeepsDefaultAlphaRegardlessOfData) {
    // The alpha computation loop starts at v = 1 by convention (node id 0 is reserved/unused),
    // so node 0 must always read back the initialized default of 0.95 in both directions.
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 4;

    const CSRMatrix matrix = build_csr_matrix(weights);

    EXPECT_FLOAT_EQ(matrix.forward_alpha(0), 0.95f);
    EXPECT_FLOAT_EQ(matrix.backward_alpha(0), 0.95f);
}

// ============================================================================
// 2. RELAXATION TOPOLOGIQUE ITERATIVE (relax_topology)
// ============================================================================
//
// Invariants / contract:
//   - Nodes with n_length[nid] == cfg::NODE_UNSEEN_32 do not exist and are never touched; their
//     output coordinate stays the sentinel.
//   - Nodes that exist but have zero forward AND zero backward degree (fully isolated) are filtered
//     out of the active set and never move away from their initial median_pos anchor.
//   - lambda_factor == 1.0 collapses `x_pred` to `anchor_pos` unconditionally (one_minus_lambda ==
//     0), so the result must equal the input median_pos exactly (up to uint32 rounding) regardless
//     of graph topology -- this is a strong, cheaply-checkable regression guard.
//   - max_iterations == 0 means the hot-loop body never executes: output must equal the rounded
//     input, and the reported iteration count is 0.
//   - An empty active set converges immediately: RMS displacement is 0 on the very first iteration
//     (0 < convergence_threshold whenever the threshold is positive), so exactly 1 iteration runs.

TEST_F(CfgFixture, LambdaOneIgnoresTopologyAndReturnsAnchorPositions) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 3;
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32, 10, 20};
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32, 100, 250};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/50, /*lambda_factor=*/1.0f);

    EXPECT_EQ(result[1], 100u);
    EXPECT_EQ(result[2], 250u);
}

TEST_F(CfgFixture, LambdaZeroMovesNodesTowardTopologicalBarycenter) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 1;
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32, 10, 10};
    // Deliberately mis-anchored so pure topology-driven relaxation must move node 2.
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32, 0, 1000};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/50, /*lambda_factor=*/0.0f);

    EXPECT_NE(result[2], 1000u);
    EXPECT_GT(iterations, 0u);
}

TEST_F(CfgFixture, ZeroMaxIterationsReturnsInputUnchanged) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 1;
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32, 10, 10};
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32, 0, 1000};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/0, /*lambda_factor=*/0.0f);

    EXPECT_EQ(result[1], 0u);
    EXPECT_EQ(result[2], 1000u);
    EXPECT_EQ(iterations, 0u);
}

TEST_F(CfgFixture, IsolatedNodeNeverMovesFromItsAnchor) {
    cfg::ARRAY_SIZE = 4;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 1; // node 3 stays fully isolated (no edges at all)
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32, 10, 10, 5};
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32, 0, 1000, 42};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/50, /*lambda_factor=*/0.0f);

    // Node 3 has zero forward and zero backward degree: filtered out of the active set entirely.
    EXPECT_EQ(result[3], 42u);
}

TEST_F(CfgFixture, NonexistentNodeKeepsSentinelInOutput) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_HAPLO = 1;

    std::unordered_map<uint64_t, uint32_t> weights;
    weights[pack_edge(1, 2)] = 1;
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32, 10, cfg::NODE_UNSEEN_32};
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32, 0, 999};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/50, /*lambda_factor=*/0.5f);

    EXPECT_EQ(result[2], cfg::NODE_UNSEEN_32);
}

TEST_F(CfgFixture, EmptyActiveSetConvergesOnFirstIteration) {
    cfg::ARRAY_SIZE = 1; // only the reserved index 0, no real nodes at all
    cfg::N_HAPLO = 1;

    const std::unordered_map<uint64_t, uint32_t> weights;
    const CSRMatrix matrix = build_csr_matrix(weights);

    const std::vector<uint32_t> n_length = {cfg::NODE_UNSEEN_32};
    const std::vector<uint32_t> median_pos = {cfg::NODE_UNSEEN_32};

    const auto [result, iterations] = relax_topology(
        matrix, n_length, median_pos,
        /*convergence_threshold=*/0.0001f, /*max_iterations=*/50, /*lambda_factor=*/0.5f);

    // Zero active nodes -> RMS displacement is 0 on the very first pass -> converges immediately.
    EXPECT_EQ(iterations, 1u);
}
