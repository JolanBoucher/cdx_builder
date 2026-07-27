//
// Created by Jolan on 2026-07-26.
//

#ifndef CDX_BUILDER_GRAPH_INDEXING_H
#define CDX_BUILDER_GRAPH_INDEXING_H
#include <vector>

// take calculate the midpoint position of each node in a linearized graph.
std::vector<uint32_t> calculate_midpoint(
    const std::vector<uint32_t>& relaxed_by_start_point,
    const std::vector<uint32_t>& n_length
    );

// for each compo, sort node according to their midpoint position, start position than nid
std::vector<uint32_t> assign_local_idx(
    const std::vector<uint16_t> &compo_by_node,
    const std::vector<uint32_t> &relaxed_start,
    const std::vector<uint32_t> &relaxed_midpoint
    );

#endif //CDX_BUILDER_GRAPH_INDEXING_H
