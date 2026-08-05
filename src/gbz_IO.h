/**
* @file gbz_IO.h
 * @brief I/O, component parsing, graph traversal, and statistics utilities for GBZ/GBWTGraph objects.
 *
 * This module provides functions to load, inspect, and process variation graphs (GBZ and GBWTGraph).
 * Key functionalities include:
 * - Extracting, cleaning, and biologically sorting component/contig names (e.g., standard chromosomes
 *   first, mitochondrial DNA, and unplaced contigs last).
 * - Assigning and building connected components via Union-Find (Disjoint-Set) data structures.
 * - Computing path statistics such as haplotype counts, edge co-occurrence weights, and node median offsets.
 *
 * @created Jolan on 2026-07-22.
 */

#ifndef GBZ_IO_H
#define GBZ_IO_H

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <gbwtgraph/gbz.h>
#include <cstddef>

/**
 * @brief Computes edge weights based on co-occurrences along paths in the forward direction.
 *
 * @param graph The GBWTGraph / GBZ graph.
 * @return std::unordered_map<uint64_t, uint32_t> Mapping of packed node pairs to their weights.
 */
std::unordered_map<uint64_t, uint32_t> compute_edge_weights(
    const gbwtgraph::GBWTGraph& graph
);

/**
 * @brief Computes and assigns connected components for all nodes in a GBWT graph.
 *
 * This function utilizes a Disjoint-Set (Union-Find) data structure with path halving
 * and union-by-rank optimizations. It traverses all nodes and their edges to merge
 * connected components, followed by a full path compression pass.
 *
 * @param gbwt A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @param parent A reference to the parent vector for the Union-Find structure
 *               (unseen nodes should be marked with cfg::NODE_UNSEEN_32).
 * @param children A reference to the rank/depth vector used for union-by-rank.
 * @return std::vector<uint32_t>& A reference to the modified parent vector containing
 *         the finalized component roots.
 */
std::vector<uint32_t>& assign_connected_components(
    const gbwtgraph::GBWTGraph& gbwt,
    std::vector<uint32_t>& parent,
    std::vector<uint16_t>& children
);

/**
 * @brief Computes and returns the dense component ID for each node in the graph.
 *
 * This function builds the Union-Find structure for the graph, compresses node roots,
 * and maps each unique sparse root to a dense, zero-indexed component identifier.
 * Unused or non-existent nodes retain the `NODE_UNSEEN_16` sentinel value.
 *
 * @param gbwt A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @param parent A reference to the Union-Find parent vector, updated in-place.
 * @param children A reference to the Union-Find children vector, updated in-place.
 * @return std::vector<uint16_t> A vector containing the dense component ID (idx) for each node index.
 */
std::vector<uint16_t> get_graph_components(
    const gbwtgraph::GBWTGraph& gbwt,
    std::vector<uint32_t>& parent,
    std::vector<uint16_t>& children
);

/**
 * @brief Binds connected component IDs to path-derived contig names, sorts them according
 *        to biological precedence, and updates `nid2compo` in-place so that CIDs match this order.
 *
 * Parses each path in the graph, associates it with a component ID, and maps component
 * identifiers to their biological sorting priority (e.g., standard chromosomes first,
 * mitochondrial DNA, and unplaced contigs last). Re-indexes node-to-component mappings accordingly.
 *
 * @param gbz Constant reference to the gbwtgraph::GBZ object.
 * @param nid2compo Reference to the vector mapping node indices to component IDs (modified in-place).
 * @param max_steps_to_check Maximum steps to validate per path (0 for unlimited).
 * @return std::vector<std::string> Sorted component names where index matches the new component ID.
 */
std::vector<std::string> bind_component_names(
    const gbwtgraph::GBZ& gbz,
    std::vector<uint16_t>& nid2compo,
    size_t max_steps_to_check = 0
);

/**
 * @brief Computes the number of unique haplotypes from the paths in a GBWT graph.
 *
 * This function iterates over all paths in the graph, extracts the sample and
 * haplotype identifiers following the PanSN naming convention (sample#haplotype#...),
 * and stores them in a vector. It handles paths that do not conform to the
 * convention with a fallback mechanism and a warning message. Finally, it sorts
 * and removes duplicates to return the exact count of unique haplotypes.
 *
 * @param graph A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @return size_t The total number of unique haplotypes found across all graph paths.
 */
size_t count_haplotypes(
    const gbwtgraph::GBWTGraph& graph
);

/**
 * @brief Computes the median offset for each node in a GBWT graph.
 *
 * This function processes the graph in three main steps:
 * 1. Pass One (Occurrence Counting): Iterates through all paths and steps to count node frequencies.
 * 2. Pass Two (CSR Table Population): Populates the offset table with actual path offsets.
 * 3. Pass Three (In-Place Median Computation): Calculates the median offset for each node.
 *
 * @param graph A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @return std::vector<uint32_t> A vector containing the computed median offset for each node index.
 */
std::vector<uint32_t> compute_nodes_median_offset(
    const gbwtgraph::GBWTGraph& graph
);

// --- Internal functions (Exposed for unit testing) ---
namespace detail {

    /**
     * @brief Robustly extracts and cleans a component or contig name from a path name.
     *
     * Parses PanSN-formatted names (e.g., Sample#Haplotype#Contig or Sample#Contig)
     * or falls back to heuristic regular expression matching for standard chromosome patterns.
     *
     * @param path_name The full input path or sequence name.
     * @return std::string The extracted component/contig identifier.
     */
    std::string extract_compo_robust(const std::string& path_name);

    /**
     * @brief Custom comparator to sort component or chromosome names logically.
     *
     * Orders chromosomes naturally (e.g., numeric chromosomes 1-22 first, followed
     * by sex chromosomes, mitochondrial DNA, alternative loci, and unplaced contigs last).
     *
     * @param a First component name string.
     * @param b Second component name string.
     * @return true if string 'a' should come before string 'b'.
     */
    bool compare_compo_names(const std::string& a, const std::string& b);

    /**
     * @brief Combines two 32-bit node IDs (src -> dst) into a single 64-bit key.
     *
     * Assumes all nids are valid 32-bit unsigned integers.
     *
     * @param src_nid The source node ID.
     * @param dst_nid The destination node ID.
     * @return uint64_t The packed 64-bit edge key.
     */
    inline uint64_t pack_node_pair(handlegraph::nid_t src_nid, handlegraph::nid_t dst_nid) {
        return (static_cast<uint64_t>(src_nid) << 32) | static_cast<uint64_t>(dst_nid);
    }
}

#endif // GBZ_IO_H