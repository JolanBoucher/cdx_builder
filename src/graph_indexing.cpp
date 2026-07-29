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
#include <vector>
#include "graph_indexing.h"

#include <numeric>

#include "constant.h"

std::vector<uint32_t> calculate_midpoint(
    const std::vector<uint32_t>& start_pos,
    const std::vector<uint32_t>& nodes_length)
{
    // Output vector preallocation
    std::vector<uint32_t> midpoint(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);

    // For each node in the start_pos vector
    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        const uint32_t len = nodes_length[nid];
        if (len == cfg::NODE_UNSEEN_32) continue;

        // The bitwise right shift '>> 1' mathematically corresponds to floor(length / 2).
        midpoint[nid] = start_pos[nid] + (len >> 1);
    }
    return midpoint;
}


std::vector<uint32_t> assign_local_idx(
    const std::vector<uint16_t> &compo_by_node,
    const std::vector<uint32_t> &relaxed_start,
    const std::vector<uint32_t> &relaxed_midpoint)
{
    // Initialize the output local index vector with the global unseen sentinel value
    std::vector<uint32_t> local_index(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);

    // Reserve a buffer to store active node identifiers
    std::vector<uint32_t> active_nodes;
    active_nodes.reserve(cfg::ARRAY_SIZE);

    // Validate node components and coordinates before pushing them into the active node buffer
    for (std::size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid)
    {
        // Skip nodes that have not been assigned a valid relaxed midpoint
        if (relaxed_midpoint[nid] == cfg::NODE_UNSEEN_32) continue;

        const uint32_t start = relaxed_start[nid];
        const uint16_t compo = compo_by_node[nid];

        // Ensure active nodes possess valid start coordinates and component bindings
        if (start == cfg::NODE_UNSEEN_32 ) {
            throw std::logic_error("Active node " + std::to_string(nid) +
                " has no relaxed start coordinate.");
        }
        if (compo == cfg::NODE_UNSEEN_16) {
            throw std::logic_error("Active node " + std::to_string(nid) +
                " is not attributed to any connected component.");
        }

        // Component count boundaries are handled externally; collect valid active node ID
        active_nodes.push_back(static_cast<uint32_t>(nid));
    }

    // Return early if no active nodes are found
    if (active_nodes.empty()) return local_index;

    // Extract raw pointers for optimal performance within the sorting lambda
    const uint16_t* compo_ptr = compo_by_node.data();
    const uint32_t* midpoint_ptr = relaxed_midpoint.data();
    const uint32_t* start_ptr = relaxed_start.data();

    // Sort active node IDs lexicographically based on spatial and topological properties.
    // The comparator evaluates criteria sequentially until a tie is broken.
    std::sort(
        active_nodes.begin(),
        active_nodes.end(),
        [compo_ptr, midpoint_ptr, start_ptr](uint32_t left, uint32_t right) noexcept {
            // Primary criterion: Group nodes by connected component ID
            const uint16_t left_compo = compo_ptr[left];
            const uint16_t right_compo = compo_ptr[right];
            if (left_compo != right_compo) {
                return left_compo < right_compo;
            }

            // Secondary criterion: Order by relaxed midpoint coordinate
            const uint32_t left_midpoint = midpoint_ptr[left];
            const uint32_t right_midpoint = midpoint_ptr[right];
            if (left_midpoint != right_midpoint) {
                return left_midpoint < right_midpoint;
            }

            // Tertiary criterion: Order by relaxed start coordinate
            const uint32_t left_start = start_ptr[left];
            const uint32_t right_start = start_ptr[right];
            if (left_start != right_start) {
                return left_start < right_start;
            }

            // Fallback criterion: Strict node ID ordering for total determinism
            return left < right;
        }
    );

    // Attribute sequential local ranks within each component group
    uint16_t curr_compo = compo_ptr[active_nodes.front()];
    uint32_t local_rank = 0;

    for (const uint32_t nid : active_nodes) {
        const uint16_t compo = compo_ptr[nid];

        // Reset local rank counter when transitioning to a new component
        if (compo != curr_compo) {
            curr_compo = compo;
            local_rank = 0;
        }

        // Guard against uint32_t index overflow for excessively large components
        if (local_rank == cfg::NODE_UNSEEN_32) {
            throw std::overflow_error("A connected component contains too many nodes "
                                      "for uint32_t local indices.");
        }

        // Assign local rank and increment counter
        local_index[nid] = local_rank;
        ++local_rank;
    }

    return local_index;
}
