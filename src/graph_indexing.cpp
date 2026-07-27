//
// Created by Jolan on 2026-07-26.
//
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

//TODO (can be parallelized if slow)
std::vector<uint32_t> assign_local_idx(
    const std::vector<uint16_t> &compo_by_node,
    const std::vector<uint32_t> &relaxed_start,
    const std::vector<uint32_t> &relaxed_midpoint)
{
   // we don't need to the vector size verification it's guaranty
   // since all vector are initialized with cfg::ARRAY_SIZE

    // we don't need to check overflow here
    // it's supposed to be guarantied by the function building those vector

    // initialize the output vector
    std::vector<uint32_t> local_index(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);

    // reserve a buffer of active node
    std::vector<uint32_t> active_nodes;
    active_nodes.reserve(cfg::ARRAY_SIZE);

    // validation of node compo and coordinate before putting in the active node buffer
    for (std::size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid)
    {
        if (relaxed_midpoint[nid] == cfg::NODE_UNSEEN_32) continue;

        const uint32_t start = relaxed_start[nid];
        const uint16_t compo = compo_by_node[nid];

        if (start == cfg::NODE_UNSEEN_32 ) {
            throw std::logic_error("Active node " + std:: to_string(nid) +
                " has no relaxed start coordinate.");
        }
        if (compo == cfg::NODE_UNSEEN_16) {
            throw std::logic_error("Active node " + std::to_string(nid) +
                " is not attributed to any connected component.");
        }

        // checking the number of component should not be checked in this function
        active_nodes.push_back(static_cast<uint32_t>(nid));
    }
    if (active_nodes.empty()) return local_index;

    // extracting the raw pointers
    const uint16_t* compo_ptr = compo_by_node.data();
    const uint32_t* midpoint_ptr = relaxed_midpoint.data();
    const uint32_t* start_ptr = relaxed_start.data();

    // Sort active node IDs lexicographically based on their spatial properties.
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

    // Attributing local rank
    uint16_t curr_compo = compo_ptr[active_nodes.front()];
    uint32_t local_rank = 0;

    for (const uint32_t nid : active_nodes) {
        const uint16_t compo = compo_ptr[nid];

        if (compo != curr_compo) {
            curr_compo = compo;
            local_rank = 0;
        }
        // Security against overflow
        if (local_rank == cfg::NODE_UNSEEN_32) {
            throw std::overflow_error("A connected component contains too many nodes"
                                      " for uint32_t local indices");
        }
        local_index[nid] = local_rank;
        ++local_rank;
    }

    return local_index;
}
