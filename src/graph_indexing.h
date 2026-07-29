/**
 * @file graph_indexing.h
 * @brief Graph Indexing & Spatial Ordering Module.
 *
 * This module provides utilities for spatial indexing and coordinate organization of graph nodes.
 * It features functions to compute node midpoint coordinates from relaxed start positions and sequence lengths,
 * and to assign contiguous local rank indices to active nodes sorted lexicographically by connected component,
 * midpoint coordinate, and start coordinate.
 */

#ifndef CDX_BUILDER_GRAPH_INDEXING_H
#define CDX_BUILDER_GRAPH_INDEXING_H
#include <vector>

/**
 * @brief Computes the midpoint coordinate for each graph node based on its start position and length.
 *
 * Iterates through all nodes up to the global array size capacity, calculating the midpoint
 * by adding half of the node's length to its start position. Nodes marked with the unseen
 * sentinel (`NODE_UNSEEN_32`) in length are skipped.
 *
 * @param relaxed_by_start_point Vector mapping each node ID to its starting coordinate position.
 * @param n_length Vector mapping each node ID to its sequence length.
 * @return std::vector<uint32_t> Vector containing calculated midpoint coordinates for all valid nodes
 *         (with unseen nodes retaining `NODE_UNSEEN_32`).
 */
std::vector<uint32_t> calculate_midpoint(
    const std::vector<uint32_t>& relaxed_by_start_point,
    const std::vector<uint32_t>& n_length
    );

/**
 * @brief Assigns local sequential rank indices to active graph nodes within their connected components.
 *
 * Filters active nodes based on valid relaxed coordinates and component assignments,
 * then sorts them lexicographically (by component ID, relaxed midpoint, relaxed start coordinate,
 * and fallback node ID). Finally, assigns a contiguous 0-based local index rank to nodes
 * grouped within each connected component.
 *
 * @param compo_by_node Vector mapping each node ID to its connected component ID.
 * @param relaxed_start Vector mapping each node ID to its relaxed start coordinate.
 * @param relaxed_midpoint Vector mapping each node ID to its relaxed midpoint coordinate.
 * @return std::vector<uint32_t> Vector mapping node IDs to their assigned local rank index
 *         (or `NODE_UNSEEN_32` for inactive/unseen nodes).
 * @throws std::logic_error If an active node lacks a valid relaxed start coordinate or component assignment.
 * @throws std::overflow_error If a connected component exceeds `uint32_t` capacity for local node indices.
 */
std::vector<uint32_t> assign_local_idx(
    const std::vector<uint16_t> &compo_by_node,
    const std::vector<uint32_t> &relaxed_start,
    const std::vector<uint32_t> &relaxed_midpoint
    );

#endif //CDX_BUILDER_GRAPH_INDEXING_H
