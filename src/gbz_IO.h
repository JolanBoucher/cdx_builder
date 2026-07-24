//
// Created by Jolan on 2026-07-22.
//

#ifndef GBZ_IO_H
#define GBZ_IO_H

#include <string>
#include <vector>
#include <gbwtgraph/gbz.h>

/**
 * @brief Charge un fichier GBZ et retourne la liste triée des composants.
 */
std::vector<std::string> get_sorted_components(
    const gbwtgraph::GBZ& gbz
    );


// Calcule la fréquence de passage dans le sens forward pour chaque arête du graphe.
std::unordered_map<uint64_t, uint32_t> compute_edge_weights(
    const gbwtgraph::GBWTGraph &graph
);

// assign each node to a connected component via union-find
std::vector<uint32_t>& assign_connected_components(
    const gbwtgraph::GBWTGraph& gbwt,
    std::vector<uint32_t>& parent,
    std::vector<uint16_t>& children
);

// assign an idx by components to each node
std::vector<uint16_t> get_graph_components(
    const gbwtgraph::GBWTGraph& gbwt,
    std::vector<uint32_t>& parent,
    std::vector<uint16_t>& children
);

// bind the name of each connected component extracted from the paths the right component id
std::vector<std::string> bind_component_names(
    const gbwtgraph::GBZ& gbz,
    const std::vector<uint16_t>& nid2compo,
    size_t max_steps_to_check = 0);

// Count and return an estimate of the number of haplotype in the graph
size_t count_haplotypes(
    const gbwtgraph::GBWTGraph& graph);


// Computes and returns the median offset for each node in the GBWT graph.
std::vector<uint32_t> compute_nodes_median_offset(
    const gbwtgraph::GBWTGraph& graph);

// --- Fonctions internes (Exposées pour les tests unitaires) ---
namespace detail {
    std::string extract_compo_robust(const std::string& path_name);

    bool compare_compo_names(const std::string& a, const std::string& b);

    /**
     * @brief Combine deux node IDs de 32 bits (src -> dst) en une seule clé de 64 bits.
     *  assume all nid are valid 32 bits uint
     */
    inline uint64_t pack_node_pair(handlegraph::nid_t src_nid, handlegraph::nid_t dst_nid) {
        return (static_cast<uint64_t>(src_nid) << 32) | static_cast<uint64_t>(dst_nid);
    }
}

#endif // GBZ_IO_H