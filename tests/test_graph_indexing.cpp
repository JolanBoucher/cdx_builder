/**
 * @file test_graph_indexing.cpp
 * @brief Unit test suite for the graph indexing module (src/graph_indexing.h/.cpp).
 *
 * Covers, in order:
 *   - calculate_midpoint (per-node midpoint coordinate from start position + sequence length)
 *   - assign_local_idx   (per-component local rank assignment via lexicographic sort)
 *
 * This module is pure and has zero external dependencies (no gbwt/gbwtgraph), which is why
 * tests/CMakeLists.txt links test_graph_indexing only against graph_indexing.cpp.
 *
 * Shared test infrastructure:
 *   - CfgFixture: saves/restores cfg::ARRAY_SIZE around each test (both functions under test size
 *     their output strictly off of it; see test_gbz_io.cpp for the full rationale of why this
 *     matters across a shared gtest binary).
 */

#include <gtest/gtest.h>
#include <cstdint>
#include <limits>
#include <vector>
#include "graph_indexing.h"
#include "constant.h"

namespace {

/// Saves/restores the cfg:: globals this module reads (cfg::ARRAY_SIZE).
class CfgFixture : public ::testing::Test {
protected:
    void SetUp() override {
        saved_array_size_ = cfg::ARRAY_SIZE;
    }

    void TearDown() override {
        cfg::ARRAY_SIZE = saved_array_size_;
    }

private:
    size_t saved_array_size_{};
};

} // namespace

// ============================================================================
// 1. CALCUL DU POINT MEDIAN PAR NOEUD (calculate_midpoint)
// ============================================================================
//
// Invariants / contract:
//   - midpoint[nid] = start_pos[nid] + floor(length[nid] / 2), via the bit-shift `len >> 1`.
//   - A node whose length is cfg::NODE_UNSEEN_32 (does not exist) is skipped and keeps the
//     cfg::NODE_UNSEEN_32 sentinel in the output, regardless of whatever start_pos holds for it.
//   - Unlike compute_alpha_factors() in graph_linearization.cpp (which starts at index 1 because
//     node id 0 is conventionally reserved/unused), this loop starts at index 0 -- node 0 IS
//     processed here if n_length[0] happens to be a valid length. Worth being explicit about,
//     since assuming the same "index 0 is always skipped" convention across modules would be wrong.
//   - No overflow protection: `start_pos + (len >> 1)` is unchecked uint32_t arithmetic. A crafted
//     start_pos close to UINT32_MAX combined with a nonzero length silently wraps around instead of
//     erroring, and could in principle land exactly on cfg::NODE_UNSEEN_32 (0xFFFFFFFF) for a node
//     that legitimately exists -- which would then be misread downstream as "does not exist".

TEST_F(CfgFixture, EvenLengthMidpointIsExactHalf) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint32_t> start_pos = {0, 100};
    const std::vector<uint32_t> length    = {cfg::NODE_UNSEEN_32, 10};

    const auto midpoint = calculate_midpoint(start_pos, length);

    EXPECT_EQ(midpoint[1], 105u); // 100 + 10/2
}

TEST_F(CfgFixture, OddLengthMidpointFloorsDown) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint32_t> start_pos = {0, 100};
    const std::vector<uint32_t> length    = {cfg::NODE_UNSEEN_32, 7};

    const auto midpoint = calculate_midpoint(start_pos, length);

    EXPECT_EQ(midpoint[1], 103u); // 100 + floor(7/2) = 100 + 3
}

TEST_F(CfgFixture, ZeroLengthMidpointEqualsStart) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint32_t> start_pos = {0, 50};
    const std::vector<uint32_t> length    = {cfg::NODE_UNSEEN_32, 0};

    const auto midpoint = calculate_midpoint(start_pos, length);

    EXPECT_EQ(midpoint[1], 50u);
}

TEST_F(CfgFixture, NonexistentNodeKeepsSentinelRegardlessOfStartPos) {
    cfg::ARRAY_SIZE = 2;
    // start_pos[1] is garbage-but-plausible; length[1] marks the node as nonexistent.
    const std::vector<uint32_t> start_pos = {0, 12345};
    const std::vector<uint32_t> length    = {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32};

    const auto midpoint = calculate_midpoint(start_pos, length);

    EXPECT_EQ(midpoint[1], cfg::NODE_UNSEEN_32);
}

TEST_F(CfgFixture, NodeIndexZeroIsProcessedUnlikeTheAlphaFactorConvention) {
    // Contrast with compute_alpha_factors() in graph_linearization.cpp, which starts at v = 1.
    // calculate_midpoint() has no such carve-out: index 0 is computed like any other index.
    cfg::ARRAY_SIZE = 1;
    const std::vector<uint32_t> start_pos = {40};
    const std::vector<uint32_t> length    = {8};

    const auto midpoint = calculate_midpoint(start_pos, length);

    EXPECT_EQ(midpoint[0], 44u); // 40 + 8/2, NOT left at the NODE_UNSEEN_32 default
}

TEST_F(CfgFixture, StartPosNearUint32MaxSilentlyWrapsInsteadOfErroring) {
    // Documents a real, unguarded overflow: no validation catches start_pos + len/2 exceeding
    // uint32_t capacity. This is well-defined (modulo 2^32) unsigned wraparound, not UB, but it is
    // silent data corruption from the caller's point of view.
    cfg::ARRAY_SIZE = 1;
    const std::vector<uint32_t> start_pos = {std::numeric_limits<uint32_t>::max() - 2};
    const std::vector<uint32_t> length    = {8}; // + 4 wraps past UINT32_MAX

    const auto midpoint = calculate_midpoint(start_pos, length);

    // (UINT32_MAX - 2) + 4, wrapped modulo 2^32, == 1.
    EXPECT_EQ(midpoint[0], 1u);
}

// ============================================================================
// 2. RANG LOCAL PAR COMPOSANTE (assign_local_idx)
// ============================================================================
//
// Invariants / contract:
//   - Only nodes with a valid (non-sentinel) relaxed_midpoint are considered "active"; everything
//     else keeps cfg::NODE_UNSEEN_32 in the output and never consumes a rank slot.
//   - An active node with an unset relaxed_start or component throws std::logic_error -- this is
//     treated as corrupt upstream data, not a recoverable case.
//   - Active nodes are sorted lexicographically by (component, midpoint, start, node id) -- the
//     node id is a pure tie-breaker for total determinism, only reached when component, midpoint,
//     AND start are all equal.
//   - Local rank restarts at 0 for every new component encountered while walking the sorted list
//     (components are contiguous after the primary sort key).
//   - The uint32_t overflow guard (a single component exceeding ~4 billion nodes) is not
//     practically exercisable in a unit test and is not covered here.

TEST_F(CfgFixture, NoActiveNodesReturnsAllSentinelWithoutThrowing) {
    cfg::ARRAY_SIZE = 3;
    const std::vector<uint16_t> compo   = {cfg::NODE_UNSEEN_16, cfg::NODE_UNSEEN_16, cfg::NODE_UNSEEN_16};
    const std::vector<uint32_t> start   = {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32};

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[0], cfg::NODE_UNSEEN_32);
    EXPECT_EQ(local_idx[1], cfg::NODE_UNSEEN_32);
    EXPECT_EQ(local_idx[2], cfg::NODE_UNSEEN_32);
}

TEST_F(CfgFixture, SingleActiveNodeGetsRankZero) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 10};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 15};

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[1], 0u);
}

TEST_F(CfgFixture, MultipleNodesSameComponentRankedByMidpointOrder) {
    cfg::ARRAY_SIZE = 4;
    // Nodes 1,2,3 all in component 0, with midpoints out of node-id order: 3 < 1 < 2.
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0, 0, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 100, 200, 50};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 100, 200, 50};

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[3], 0u); // midpoint 50 -> first
    EXPECT_EQ(local_idx[1], 1u); // midpoint 100 -> second
    EXPECT_EQ(local_idx[2], 2u); // midpoint 200 -> third
}

TEST_F(CfgFixture, TiedMidpointBrokenByStart) {
    cfg::ARRAY_SIZE = 3;
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 200, 100}; // node 2 starts earlier
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 500, 500}; // tied midpoint

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[2], 0u); // lower start wins the tie
    EXPECT_EQ(local_idx[1], 1u);
}

TEST_F(CfgFixture, TiedMidpointAndStartBrokenByNodeId) {
    cfg::ARRAY_SIZE = 3;
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 100, 100};    // tied
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 500, 500};    // tied

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[1], 0u); // lower node id wins the final tie-break
    EXPECT_EQ(local_idx[2], 1u);
}

TEST_F(CfgFixture, MultipleComponentsRestartRankAtZeroEach) {
    cfg::ARRAY_SIZE = 5;
    // Nodes 1,2 in component 0; nodes 3,4 in component 1.
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0, 0, 1, 1};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 10, 20, 10, 20};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 10, 20, 10, 20};

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[1], 0u);
    EXPECT_EQ(local_idx[2], 1u);
    EXPECT_EQ(local_idx[3], 0u); // component 1's rank restarts at 0
    EXPECT_EQ(local_idx[4], 1u);
}

TEST_F(CfgFixture, MixedActiveAndInactiveNodesLeaveInactiveAsSentinel) {
    cfg::ARRAY_SIZE = 3;
    // Node 2 is inactive (no midpoint) even though it has a component and start assigned.
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 10, 20};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 10, cfg::NODE_UNSEEN_32};

    const auto local_idx = assign_local_idx(compo, start, midpoint);

    EXPECT_EQ(local_idx[1], 0u);              // the only active node gets rank 0
    EXPECT_EQ(local_idx[2], cfg::NODE_UNSEEN_32); // inactive: untouched, does not consume a rank
}

TEST_F(CfgFixture, ActiveNodeWithMissingStartThrowsLogicError) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32}; // missing!
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 15}; // but midpoint says "active"

    EXPECT_THROW(assign_local_idx(compo, start, midpoint), std::logic_error);
}

TEST_F(CfgFixture, ActiveNodeWithMissingComponentThrowsLogicError) {
    cfg::ARRAY_SIZE = 2;
    const std::vector<uint16_t> compo    = {cfg::NODE_UNSEEN_16, cfg::NODE_UNSEEN_16}; // missing!
    const std::vector<uint32_t> start    = {cfg::NODE_UNSEEN_32, 10};
    const std::vector<uint32_t> midpoint = {cfg::NODE_UNSEEN_32, 15}; // but midpoint says "active"

    EXPECT_THROW(assign_local_idx(compo, start, midpoint), std::logic_error);
}
