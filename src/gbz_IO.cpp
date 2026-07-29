#include "gbz_IO.h"

#include <vector>
#include <string>
#include <unordered_set>
#include <algorithm>
#include <regex>
#include <boost/xpressive/regex_primitives.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwt/metadata.h>

#include "constant.h"

namespace detail {

    static std::string clean_subrange(const std::string& name) {
        size_t bracket_pos = name.find('[');
        if (bracket_pos != std::string::npos) {
            return name.substr(0, bracket_pos);
        }
        size_t colon_pos = name.find(':');
        if (colon_pos != std::string::npos) {
            return name.substr(0, colon_pos);
        }
        return name;
    }

    std::string extract_compo_robust(const std::string& name) {
        // Parsing PanSN (Sample#Haplotype#Contig[#Index] OU Sample#Contig)
        size_t first = name.find('#');
        if (first != std::string::npos) {
            size_t second = name.find('#', first + 1);
            if (second != std::string::npos) {
                // Format à 3+ champs : Sample#Haplotype#Contig...
                size_t third = name.find('#', second + 1);
                std::string compo = (third == std::string::npos)
                                    ? name.substr(second + 1)
                                    : name.substr(second + 1, third - second - 1);
                return clean_subrange(compo);
            } else {
                // Format à 2 champs incomplet : Sample#Contig
                std::string compo = name.substr(first + 1);
                return clean_subrange(compo);
            }
        }

        // Nettoyage des sous-intervalles
        std::string cleaned_name = clean_subrange(name);

        // Recherche heuristique Regex pour les motifs de chromosomes standards
        static const std::regex chr_regex("(chr[0-9A-Za-z]+|[0-9]{1,2}|chrMT|chrM|chrmt)", std::regex_constants::icase);
        std::smatch match;
        if (std::regex_search(cleaned_name, match, chr_regex)) {
            return match.str(1);
        }

        return cleaned_name;
    }

    bool compare_compo_names(const std::string& name_a, const std::string& name_b) {
    auto is_mitochondrial = [](const std::string& s) {
        return (s.find("chrmt") != std::string::npos ||
                s.find("chrMT") != std::string::npos ||
                s.find("chrM")   != std::string::npos ||
                s.find("chrm")   != std::string::npos ||
                s == "MT" || s == "M");
    };

    auto is_unplaced = [](const std::string& s) {
        return (s.find("chrUn") != std::string::npos ||
                s.find("random")!= std::string::npos ||
                s.find("dec")   != std::string::npos ||
                s.find("Un_")   != std::string::npos);
    };

    auto is_alt = [](const std::string& s) {
        return (s.find("alt") != std::string::npos);
    };

    // 1. Les contigs non placés (chrUn) vont TOUT à la fin (après chrMT)
    bool is_un_a = is_unplaced(name_a);
    bool is_un_b = is_unplaced(name_b);
    if (is_un_a != is_un_b) {
        return is_un_b; // Si 'b' est unplaced, 'a' (ex: chrMT) passe devant !
    }

    // 2. La mitochondrie vient juste avant les unplaced (mais après les chromosomes).
    bool is_mt_a = is_mitochondrial(name_a);
    bool is_mt_b = is_mitochondrial(name_b);
    if (is_mt_a != is_mt_b) {
        return is_mt_b; // Le chromosome standard/sexuel passe avant MT
    }

    // 3. Extraction du numéro de chromosome principal (ex: "chr1", "chr1_alt" -> 1)
    auto extract_num = [](const std::string& s) {
        static const std::regex num_regex("(?:chr)?([0-9]{1,2})", std::regex_constants::icase);
        std::smatch match;
        if (std::regex_search(s, match, num_regex)) {
            return std::stoi(match.str(1));
        }
        return -1;
    };

    int num_a = extract_num(name_a);
    int num_b = extract_num(name_b);

    // Si les deux ont des numéros différents (ex : chr1_alt [1] vs chr2 [2]).
    if (num_a != num_b && num_a != -1 && num_b != -1) {
        return num_a < num_b;
    }

    // 4. Si même numéro (ex: chr1 vs chr1_alt), le canonique passe avant l'alt
    bool is_alt_a = is_alt(name_a);
    bool is_alt_b = is_alt(name_b);
    if (is_alt_a != is_alt_b) {
        return is_alt_b;
    }

    // 5. Tri par défaut
    return name_a < name_b;
}
} // namespace detail

/**
 * @brief Extraction sur l'objet GBZ avec fallback métadonnées GBWT.
 */
std::string extract_compo_name(const gbwtgraph::GBZ& gbz, const handlegraph::path_handle_t& ph) {
    const gbwt::size_type path_id = gbz.graph.handle_to_path(ph);

    if (gbz.index.hasMetadata()) {
        const gbwt::Metadata& metadata = gbz.index.metadata;

        if (metadata.hasPathNames() && metadata.hasContigNames() && path_id < metadata.paths()) {
            const gbwt::PathName& path_info = metadata.path(path_id);

            if (path_info.contig < metadata.contigs()) {
                std::string contig_name = metadata.contig(path_info.contig);
                if (!contig_name.empty()) {
                    return detail::clean_subrange(contig_name);
                }
            }
        }
    }

    std::string name = gbz.graph.get_path_name(ph);
    return detail::extract_compo_robust(name);
}


/**
 * @brief Binds connected component IDs to path-derived contig names, sorts them according
 *        to biological precedence, and updates `nid2compo` in-place so that CIDs match this order.
 *
 * @param gbz A constant reference to the gbwtgraph::GBZ object.
 * @param nid2compo Reference to the vector mapping node indices to component IDs (modified in-place).
 * @param max_steps_to_check Maximum steps to validate per path (0 for unlimited).
 * @return std::vector<std::string> Sorted component names where index == new component_id.
 */
std::vector<std::string> bind_component_names(
    const gbwtgraph::GBZ& gbz,
    std::vector<uint16_t>& nid2compo,
    const size_t max_steps_to_check)
{
    std::vector<std::string> raw_component_names(cfg::N_COMPO);

    // 1. Assign raw component names based on GBZ paths
    gbz.graph.for_each_path_handle(
        [&](const handlegraph::path_handle_t& path) -> bool {
            const std::string path_name = gbz.graph.get_path_name(path);
            const std::string component_name = extract_compo_name(gbz, path);

            uint16_t path_component = cfg::NODE_UNSEEN_16;
            size_t steps_checked = 0;

            gbz.graph.for_each_step_in_path(
                path,
                [&](const handlegraph::step_handle_t& step) -> bool {
                    if (max_steps_to_check > 0 && steps_checked >= max_steps_to_check) {
                        return false;
                    }

                    const handlegraph::handle_t handle = gbz.graph.get_handle_of_step(step);
                    const handlegraph::nid_t nid = gbz.graph.get_id(handle);
                    const auto node_index = static_cast<size_t>(nid);

                    if (node_index >= nid2compo.size()) {
                        throw std::runtime_error("Node ID outside nid2compo: " + std::to_string(nid));
                    }

                    const uint16_t component_id = nid2compo[node_index];

                    if (component_id == cfg::NODE_UNSEEN_16) {
                        throw std::runtime_error("Path " + path_name + " contains a node without a component");
                    }

                    if (component_id >= raw_component_names.size()) {
                        throw std::runtime_error("Invalid component ID " + std::to_string(component_id) +
                                                 " for node " + std::to_string(nid));
                    }

                    if (path_component == cfg::NODE_UNSEEN_16) {
                        path_component = component_id;
                    }
                    else if (path_component != component_id) {
                        throw std::runtime_error("Path " + path_name + " crosses multiple connected components");
                    }

                    steps_checked++;
                    return true;
                }
            );

            if (path_component == cfg::NODE_UNSEEN_16) {
                return true;
            }

            std::string& assigned_name = raw_component_names[path_component];

            if (assigned_name.empty()) {
                assigned_name = component_name;
            }
            else if (assigned_name != component_name) {
                throw std::runtime_error(
                    "Connected component " + std::to_string(path_component) +
                    " has conflicting names: " + assigned_name + " and " + component_name
                );
            }
            return true;
        }
    );

    // 2. Validate that every connected component was assigned a path name
    for (size_t component_id = 0; component_id < raw_component_names.size(); ++component_id) {
        if (raw_component_names[component_id].empty()) {
            throw std::runtime_error(
                "Connected component " + std::to_string(component_id) +
                " has no associated path name"
            );
        }
    }

    // 3. Pre-extract/normalize component names to avoid calling std::regex inside std::sort
    std::vector<std::string> normalized_names(cfg::N_COMPO);
    for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
        normalized_names[cid] = detail::extract_compo_robust(raw_component_names[cid]);
    }

    // 4. Compute permutation vector based on biological sorting rules
    std::vector<uint16_t> order(cfg::N_COMPO);
    std::iota(order.begin(), order.end(), static_cast<uint16_t>(0));

    std::sort(order.begin(), order.end(), [&](uint16_t a, uint16_t b) {
        return detail::compare_compo_names(normalized_names[a], normalized_names[b]);
    });

    // 5. Build inverse lookup table (old_cid -> new_cid)
    std::vector<uint16_t> cid_remap(cfg::N_COMPO);
    for (uint16_t new_cid = 0; new_cid < cfg::N_COMPO; ++new_cid) {
        cid_remap[order[new_cid]] = new_cid;
    }

    // 6. In-place update of nid2compo to match biological CIDs
    for (uint16_t& cid : nid2compo) {
        if (cid == cfg::NODE_UNSEEN_16) continue;
        cid = cid_remap[cid];
    }

    // 7. Reorder raw component names according to the new CIDs
    std::vector<std::string> sorted_component_names(cfg::N_COMPO);
    for (uint16_t old_cid = 0; old_cid < cfg::N_COMPO; ++old_cid) {
        sorted_component_names[cid_remap[old_cid]] = std::move(raw_component_names[old_cid]);
    }

#ifndef NDEBUG
    // Optional debug printing to stdout/stderr to verify mapping in dev builds
    for (size_t cid = 0; cid < sorted_component_names.size(); ++cid) {
        std::cerr << "[DEBUG] CID " << cid << " -> " << sorted_component_names[cid] << '\n';
    }
#endif

    return sorted_component_names;
}


/**
 * @brief Prend l'objet GBZ et renvoie un vecteur trié de noms de composants uniques.
 */
std::vector<std::string> get_sorted_components(const gbwtgraph::GBZ& gbz) {
    std::unordered_set<std::string> unique_components;

    gbz.graph.for_each_path_handle([&](const handlegraph::path_handle_t& ph) {
        unique_components.insert(extract_compo_name(gbz, ph));
    });

    std::vector<std::string> components(unique_components.begin(), unique_components.end());
    std::sort(components.begin(), components.end(), detail::compare_compo_names);

    return components;
}


/**
 * @brief Computes and assigns connected components for all nodes in a GBWT graph.
 *
 * This function utilizes a Disjoint-Set (Union-Find) data structure with path halving
 * and union-by-rank optimizations. It traverses all nodes and their edges to merge
 * connected components, followed by a full path compression pass.
 *
 * @param graph A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @param parent A reference to the parent vector for the Union-Find structure
 *               (unseen nodes should be marked with cfg::NODE_UNSEEN_32).
 * @param children A reference to the rank/depth vector used for union-by-rank.
 * @return std::vector<uint32_t>& A reference to the modified parent vector containing
 *         the finalized component roots.
 */
std::vector<uint32_t> & assign_connected_components(
    const gbwtgraph::GBWTGraph &graph,
    std::vector<uint32_t> &parent,
    std::vector<uint16_t> &children)  {

    auto union_nodes = [&](handlegraph::nid_t source, const handlegraph::nid_t destination){
        auto a = static_cast<uint32_t>(source);

        // Find root of source with path halving
        while (parent[a] != a) {
            parent[a] = parent[parent[a]];
            a = parent[a];
        }
        const uint32_t root_a = a;

        auto b = static_cast<uint32_t>(destination);
        // Find root of destination with path halving
        while (parent[b] != b) {
            parent[b] = parent[parent[b]];
            b = parent[b];
        }
        const uint32_t root_b = b;

        // Union by rank
        if (root_a != root_b) {
            if (children[root_a] < children[root_b]) {
                parent[root_a] = root_b;
            }
            else if (children[root_a] > children[root_b]) {
                parent[root_b] = root_a;
            }
            else {
                parent[root_b] = root_a;
                children[root_a]++;
            }
        }
    };

    graph.for_each_handle([&](const handlegraph::handle_t& source_handle) -> bool {
        const handlegraph::nid_t source = graph.get_id(source_handle);

        // follow edges on the right side
        graph.follow_edges(source_handle, false,
            [&](const handlegraph::handle_t& destination_handle) -> bool {
                const handlegraph::nid_t destination = graph.get_id(destination_handle);
                union_nodes(source, destination);
                return true;
        });

        // Follow edges on the left side
        graph.follow_edges(source_handle, true,
            [&](const handlegraph::handle_t& destination_handle) -> bool {
                const handlegraph::nid_t destination = graph.get_id(destination_handle);
                union_nodes(source, destination);
                return true;
        });

        return true; // needed by for_each_handle
    });

    // Final compression
    for (size_t nid = 0; nid < parent.size(); ++nid) {
        if (parent[nid] == cfg::NODE_UNSEEN_32) {
            continue;
        }
        auto node = static_cast<uint32_t>(nid);
        while (parent[node] != node) {
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        parent[nid] = node;
    }

    return parent;
}


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
    std::vector<uint16_t>& children)
{
    // Build the Union-Find structure and compress the roots.
    assign_connected_components(gbwt, parent, children);

    std::vector<uint16_t> compo(cfg::ARRAY_SIZE,cfg::NODE_UNSEEN_16);
    std::unordered_map<uint32_t, uint16_t> root_to_component;

    root_to_component.reserve(1024);     // The expected number of components is small.
    uint32_t n_compo = 0;

    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        const uint32_t root = parent[nid];

        if (root == cfg::NODE_UNSEEN_32) continue; // Skip node IDs that do not exist in the graph.
        auto component = root_to_component.find(root);

        if (component == root_to_component.end()) {
            if (n_compo >= cfg::NODE_UNSEEN_16) { // NODE_UNSEEN_16 is reserved and cannot be a valid component ID.
                throw std::overflow_error("The number of connected components exceeds uint16_t capacity");
            }

            const auto component_id = static_cast<uint16_t>(n_compo);
            auto [fst, snd] = root_to_component.emplace(root, component_id);
            component = fst;
            n_compo++;
        }

        compo[nid] = component->second;
    }

    cfg::N_COMPO = n_compo;
    return compo;
}

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
size_t count_haplotypes(const gbwtgraph::GBWTGraph& graph)
{
    std::vector<std::pair<std::string, std::string>> haplotypes;
    haplotypes.reserve(graph.get_path_count());

    graph.for_each_path_handle([&](const handlegraph::path_handle_t& path) -> bool {
        std::string path_name = graph.get_path_name(path);

        const size_t first_separator = path_name.find('#');

        if (first_separator == std::string::npos) {
            haplotypes.emplace_back(path_name, "0");

            std::cerr
                << "Warning: Path name "
                << path_name
                << " does not follow the PanSN convention. "
                << "Path considered as one haplotype.\n";

            return true;
        }

        const size_t second_separator = path_name.find(
            '#',
            first_separator + 1
        );

        std::string sample = path_name.substr(
            0,
            first_separator
        );

        std::string haplotype = path_name.substr(
            first_separator + 1,
            second_separator == std::string::npos
                ? std::string::npos
                : second_separator - first_separator - 1
        );

        haplotypes.emplace_back(
            std::move(sample),
            std::move(haplotype)
        );

        return true;
    });

    std::sort(
        haplotypes.begin(),
        haplotypes.end()
    );

    const auto unique_end = std::unique(
        haplotypes.begin(),
        haplotypes.end()
    );

    return static_cast<size_t>(
        std::distance(
            haplotypes.begin(),
            unique_end
        )
    );
}


/**
 * @brief Computes the median offset for each node in a GBWT graph.
 *
 * This function processes the graph in three main steps:
 *
 * 1. **Pass One (Occurrence Counting):** Iterates through all paths and steps
 *    in the graph to count the frequency of each node, then transforms these
 *    counts into cumulative CSR (Compressed Sparse Row) offsets (`offsets_key`).
 *
 * 2. **Pass Two (CSR Table Population):** Re-iterates through the paths to populate
 *    the main `offsets_table` with the actual path offsets for each node occurrence,
 *    using an active cursor tracker.
 *
 * 3. **Pass Three (In-Place Median Computation):** Calculates the median offset
 *    for each node slice in the table. It handles edge cases (slices of size 1 or 2)
 *    and uses partial sorting (`std::nth_element`) to efficiently determine the
 *    median value in linear time for larger slices.
 *
 * @param graph A constant reference to the gbwtgraph::GBWTGraph to analyze.
 * @return std::vector<uint32_t> A vector containing the computed median offset for each node index.
 */
std::vector<uint32_t> compute_nodes_median_offset(const gbwtgraph::GBWTGraph& graph)
{
    std::vector<uint64_t> offsets_key(cfg::ARRAY_SIZE + 1, 0);

    // --- 1. Pass One: Count node occurrences ---
    graph.for_each_path_handle([&](const handlegraph::path_handle_t& path) -> bool {
        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t& step) -> bool {
            handlegraph::handle_t handle = graph.get_handle_of_step(step);
            handlegraph::nid_t nid = graph.get_id(handle);

            offsets_key[static_cast<size_t>(nid) + 1]++;
            return true;
        });
        return true;
    });

    // Transform occurrence counts into CSR offsets.
    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        offsets_key[nid + 1] += offsets_key[nid];
    }

    const auto total_offsets = static_cast<size_t>(offsets_key.back());
    std::vector<uint32_t> offsets_table(total_offsets,cfg::NODE_UNSEEN_32);
    std::vector<uint64_t> offsets_cursor(offsets_key.begin(), offsets_key.end() - 1);

    // --- 2. Pass Two: Fill the CSR offset table ---
    graph.for_each_path_handle([&](const handlegraph::path_handle_t& path) -> bool {
        uint64_t path_offset = 0;

        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t& step) -> bool {
            handlegraph::handle_t handle = graph.get_handle_of_step(step);
            handlegraph::nid_t nid = graph.get_id(handle);

            const auto node_index = static_cast<size_t>(nid);
            const auto table_index = static_cast<size_t>(offsets_cursor[node_index]++);

            if (path_offset >= cfg::NODE_UNSEEN_32) {
                throw std::overflow_error("Path offset exceeds uint32_t capacity");
            }
            offsets_table[table_index] = static_cast<uint32_t>(path_offset);
            path_offset += static_cast<uint64_t>(graph.get_length(handle));
            return true;
        });

        return true;
    });

    // offsets_cursor is no longer required.
    offsets_cursor.clear();
    offsets_cursor.shrink_to_fit();

    // --- 3. Compute medians in place ---
    // A nonexistent node or a node without any path occurrence keeps the value NODE_UNSEEN_32.
    std::vector<uint32_t> medians(cfg::ARRAY_SIZE,cfg::NODE_UNSEEN_32);

    for (size_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        const auto start = static_cast<size_t>(offsets_key[nid]);
        const auto end = static_cast<size_t>(offsets_key[nid + 1]);
        const size_t slice_size = end - start;

        // An absent node or a node without a path occurrence
        // keeps NODE_UNSEEN_32.
        if (slice_size == 0) continue;

        auto slice_begin =offsets_table.begin() + start;
        const auto slice_end = offsets_table.begin() + end;
        auto right_median = slice_begin + (slice_size / 2);

        /*
         * Select the same value as:  sorted_offsets[slice_size / 2]
         * This is the right median for an even number of observations.
         */
        std::nth_element(slice_begin, right_median,slice_end);
        medians[nid] = *right_median;
    }

    return medians;
}

//TODO Chek if their is no handle that can do that instead of this weird way
/**
 * @brief Computes edge weights based on co-occurrences along paths.
 *
 * @param graph The GBWTGraph / GBZ graph.
 * @return std::unordered_map<uint64_t, uint32_t> Mapping of packed node pairs to their weights.
 */
std::unordered_map<uint64_t, uint32_t> compute_edge_weights(const gbwtgraph::GBWTGraph& graph)
{
    std::unordered_map<uint64_t, uint32_t> weights;
    weights.reserve(graph.get_edge_count());

    graph.for_each_path_handle([&](const handlegraph::path_handle_t& path) -> bool {
        handlegraph::nid_t previous_nid = 0;
        bool has_previous_node = false;

        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t& step) -> bool {
            handlegraph::handle_t handle = graph.get_handle_of_step(step);
            handlegraph::nid_t nid = graph.get_id(handle);

            if (has_previous_node) {
                uint64_t edge_key = detail::pack_node_pair(previous_nid, nid);
                weights[edge_key]++;
            }

            previous_nid = nid;
            has_previous_node = true;
            return true;
        });
        return true;
    });

    return weights;
}