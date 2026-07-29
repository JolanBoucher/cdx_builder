#include "gbz_IO.h"
#include <vector>
#include <string>
#include <algorithm>
#include <regex>
#include <boost/xpressive/regex_primitives.hpp>
#include <gbwtgraph/gbz.h>
#include <gbwt/metadata.h>
#include "constant.h"

namespace detail {
    /**
     * @brief Cleans a sequence or path name by removing any subrange specifiers.
     *
     * Strips out bracket notations (e.g., "[...]") and colon-based coordinate
     * ranges (e.g., ":start-end") from the input string to extract the base name.
     *
     * @param name The original path or contig name string.
     * @return std::string The base name without subrange suffixes.
     */
    static std::string clean_subrange(const std::string &name) {
        // Strip trailing bracket subrange specifiers (e.g., chr1[100-200] -> chr1)
        const size_t bracket_pos = name.find('[');
        if (bracket_pos != std::string::npos) {
            return name.substr(0, bracket_pos);
        }

        // Strip trailing colon coordinate ranges (e.g., chr1:100-200 -> chr1)
        const size_t colon_pos = name.find(':');
        if (colon_pos != std::string::npos) {
            return name.substr(0, colon_pos);
        }

        return name;
    }

    std::string extract_compo_robust(const std::string &path_name) {
        // Parse PanSN naming format (Sample#Haplotype#Contig[#Index] or Sample#Contig)
        const size_t first = path_name.find('#');
        if (first != std::string::npos) {
            const size_t second = path_name.find('#', first + 1);
            if (second != std::string::npos) {
                // Format with 3+ fields: Sample#Haplotype#Contig...
                const size_t third = path_name.find('#', second + 1);
                const std::string compo = third == std::string::npos
                                              ? path_name.substr(second + 1)
                                              : path_name.substr(second + 1, third - second - 1);
                return clean_subrange(compo);
            }
            // Incomplete 2-field format: Sample#Contig
            const std::string compo = path_name.substr(first + 1);
            return clean_subrange(compo);
        }

        // Clean any subrange suffixes from the input path name
        std::string cleaned_name = clean_subrange(path_name);

        // Heuristic regex search to match standard chromosome patterns
        static const std::regex chr_regex("(chr[0-9A-Za-z]+|[0-9]{1,2}|chrMT|chrM|chrmt)", std::regex_constants::icase);
        std::smatch match;
        if (std::regex_search(cleaned_name, match, chr_regex)) {
            return match.str(1);
        }

        return cleaned_name;
    }


    bool compare_compo_names(const std::string &a, const std::string &b) {
        auto is_mitochondrial = [](const std::string &s) {
            return s.find("chrmt") != std::string::npos ||
                   s.find("chrMT") != std::string::npos ||
                   s.find("chrM") != std::string::npos ||
                   s.find("chrm") != std::string::npos ||
                   s == "MT" || s == "M";
        };

        auto is_unplaced = [](const std::string &s) {
            return s.find("chrUn") != std::string::npos ||
                   s.find("random") != std::string::npos ||
                   s.find("dec") != std::string::npos ||
                   s.find("Un_") != std::string::npos;
        };

        auto is_alt = [](const std::string &s) {
            return s.find("alt") != std::string::npos;
        };

        // 1. Unplaced contigs (chrUn) go entirely to the end (after chrMT)
        const bool is_un_a = is_unplaced(a);
        const bool is_un_b = is_unplaced(b);
        if (is_un_a != is_un_b) {
            return is_un_b; // If 'b' is unplaced, 'a' (e.g., chrMT) comes first
        }

        // 2. Mitochondrial DNA comes just before unplaced (but after standard chromosomes)
        const bool is_mt_a = is_mitochondrial(a);
        const bool is_mt_b = is_mitochondrial(b);
        if (is_mt_a != is_mt_b) {
            return is_mt_b; // Standard/sex chromosomes come before MT
        }

        // 3. Extract the primary chromosome number (e.g., "chr1", "chr1_alt" -> 1)
        auto extract_num = [](const std::string &s) {
            static const std::regex num_regex("(?:chr)?([0-9]{1,2})", std::regex_constants::icase);
            std::smatch match;
            if (std::regex_search(s, match, num_regex)) {
                return std::stoi(match.str(1));
            }
            return -1;
        };

        const int num_a = extract_num(a);
        const int num_b = extract_num(b);

        // If both have different valid numbers (e.g., chr1 [1] vs chr2 [2])
        if (num_a != num_b && num_a != -1 && num_b != -1) {
            return num_a < num_b;
        }

        // 4. If same number (e.g., chr1 vs chr1_alt), canonical comes before alt
        const bool is_alt_a = is_alt(a);
        const bool is_alt_b = is_alt(b);
        if (is_alt_a != is_alt_b) {
            return is_alt_b;
        }

        // 5. Fallback to standard lexicographical comparison
        return a < b;
    }
} // namespace detail

/**
 * @brief Extracts a component name from a GBZ object using GBWT metadata fallback.
 *
 * Attempts to retrieve and clean the contig name directly from GBWT metadata paths.
 * If metadata is missing or incomplete, falls back to parsing the path name
 * using robust extraction heuristics.
 *
 * @constref gbz The GBZ graph and index containing metadata and paths.
 * @constref ph The path handle for the target path in the graph.
 * @return std::string The cleaned component or contig name.
 */
std::string extract_compo_name(const gbwtgraph::GBZ &gbz, const handlegraph::path_handle_t &ph) {
    const gbwt::size_type path_id = gbz.graph.handle_to_path(ph);

    // Try extracting the contig name directly from GBWT metadata if available
    if (gbz.index.hasMetadata()) {
        const gbwt::Metadata &metadata = gbz.index.metadata;

        if (metadata.hasPathNames() && metadata.hasContigNames() && path_id < metadata.paths()) {
            const gbwt::PathName &path_info = metadata.path(path_id);

            if (path_info.contig < metadata.contigs()) {
                const std::string contig_name = metadata.contig(path_info.contig);
                if (!contig_name.empty()) {
                    return detail::clean_subrange(contig_name);
                }
            }
        }
    }

    // Fallback: retrieve the full path name from the graph and extract robustly
    const std::string name = gbz.graph.get_path_name(ph);
    return detail::extract_compo_robust(name);
}


std::vector<std::string> bind_component_names(
    const gbwtgraph::GBZ &gbz,
    std::vector<uint16_t> &nid2compo,
    const size_t max_steps_to_check) {
    std::vector<std::string> raw_component_names(cfg::N_COMPO);

    // 1. Assign raw component names based on GBZ paths
    gbz.graph.for_each_path_handle(
        [&](const handlegraph::path_handle_t &path) -> bool {
            const std::string path_name = gbz.graph.get_path_name(path);
            const std::string component_name = extract_compo_name(gbz, path);

            uint16_t path_component = cfg::NODE_UNSEEN_16;
            size_t steps_checked = 0;

            gbz.graph.for_each_step_in_path(
                path,
                [&](const handlegraph::step_handle_t &step) -> bool {
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
                    } else if (path_component != component_id) {
                        throw std::runtime_error("Path " + path_name + " crosses multiple connected components");
                    }

                    steps_checked++;
                    return true;
                }
            );

            if (path_component == cfg::NODE_UNSEEN_16) {
                return true;
            }

            std::string &assigned_name = raw_component_names[path_component];

            if (assigned_name.empty()) {
                assigned_name = component_name;
            } else if (assigned_name != component_name) {
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

    std::sort(order.begin(), order.end(), [&](const uint16_t a, const uint16_t b) {
        return detail::compare_compo_names(normalized_names[a], normalized_names[b]);
    });

    // 5. Build inverse lookup table (old_cid -> new_cid)
    std::vector<uint16_t> cid_remap(cfg::N_COMPO);
    for (uint16_t new_cid = 0; new_cid < cfg::N_COMPO; ++new_cid) {
        cid_remap[order[new_cid]] = new_cid;
    }

    // 6. In-place update of nid2compo to match biological CIDs
    for (uint16_t &cid: nid2compo) {
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


std::vector<uint32_t> &assign_connected_components(
    const gbwtgraph::GBWTGraph &gbwt,
    std::vector<uint32_t> &parent,
    std::vector<uint16_t> &children) {
    auto union_nodes = [&](const handlegraph::nid_t source, const handlegraph::nid_t destination) {
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
            } else if (children[root_a] > children[root_b]) {
                parent[root_b] = root_a;
            } else {
                parent[root_b] = root_a;
                children[root_a]++;
            }
        }
    };

    gbwt.for_each_handle([&](const handlegraph::handle_t &source_handle) -> bool {
        const handlegraph::nid_t source = gbwt.get_id(source_handle);

        // follow edges on the right side
        gbwt.follow_edges(source_handle, false,
                          [&](const handlegraph::handle_t &destination_handle) -> bool {
                              const handlegraph::nid_t destination = gbwt.get_id(destination_handle);
                              union_nodes(source, destination);
                              return true;
                          });

        // Follow edges on the left side
        gbwt.follow_edges(source_handle, true,
                          [&](const handlegraph::handle_t &destination_handle) -> bool {
                              const handlegraph::nid_t destination = gbwt.get_id(destination_handle);
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


std::vector<uint16_t> get_graph_components(
    const gbwtgraph::GBWTGraph &gbwt,
    std::vector<uint32_t> &parent,
    std::vector<uint16_t> &children) {
    // Build the Union-Find structure and compress the roots.
    assign_connected_components(gbwt, parent, children);

    std::vector compo(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);
    std::unordered_map<uint32_t, uint16_t> root_to_component;

    root_to_component.reserve(1024); // The expected number of components is small.
    uint32_t n_compo = 0;

    for (size_t nid = 0; nid < cfg::ARRAY_SIZE; ++nid) {
        const uint32_t root = parent[nid];

        if (root == cfg::NODE_UNSEEN_32) continue; // Skip node IDs that do not exist in the graph.
        auto component = root_to_component.find(root);

        if (component == root_to_component.end()) {
            if (n_compo >= cfg::NODE_UNSEEN_16) {
                // NODE_UNSEEN_16 is reserved and cannot be a valid component ID.
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


size_t count_haplotypes(const gbwtgraph::GBWTGraph &graph) {
    std::vector<std::pair<std::string, std::string> > haplotypes;
    haplotypes.reserve(graph.get_path_count());

    graph.for_each_path_handle([&](const handlegraph::path_handle_t &path) -> bool {
        std::string path_name = graph.get_path_name(path);

        const size_t first_separator = path_name.find('#');

        if (first_separator == std::string::npos) {
            haplotypes.emplace_back(path_name, "0");

            std::cerr << "Warning: Path name " << path_name << " does not follow the PanSN convention. "
                    << "Path considered as one haplotype.\n";

            return true;
        }

        const size_t second_separator = path_name.find('#', first_separator + 1);

        std::string sample = path_name.substr(0, first_separator);
        std::string haplotype = path_name.substr(first_separator + 1,
                                                 second_separator == std::string::npos
                                                     ? std::string::npos
                                                     : second_separator - first_separator - 1);

        haplotypes.emplace_back(std::move(sample), std::move(haplotype));
        return true;
    });

    std::sort(haplotypes.begin(), haplotypes.end());
    const auto unique_end = std::unique(haplotypes.begin(), haplotypes.end());
    return static_cast<size_t>(std::distance(haplotypes.begin(), unique_end));
}


std::vector<uint32_t> compute_nodes_median_offset(const gbwtgraph::GBWTGraph &graph) {
    std::vector<uint64_t> offsets_key(cfg::ARRAY_SIZE + 1, 0);

    // --- 1. Pass One: Count node occurrences ---
    graph.for_each_path_handle([&](const handlegraph::path_handle_t &path) -> bool {
        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t &step) -> bool {
            const handlegraph::handle_t handle = graph.get_handle_of_step(step);
            const handlegraph::nid_t nid = graph.get_id(handle);

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
    std::vector offsets_table(total_offsets, cfg::NODE_UNSEEN_32);
    std::vector offsets_cursor(offsets_key.begin(), offsets_key.end() - 1);

    // --- 2. Pass Two: Fill the CSR offset table ---
    graph.for_each_path_handle([&](const handlegraph::path_handle_t &path) -> bool {
        uint64_t path_offset = 0;

        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t &step) -> bool {
            const handlegraph::handle_t handle = graph.get_handle_of_step(step);
            const handlegraph::nid_t nid = graph.get_id(handle);

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
    std::vector medians(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);

    for (size_t nid = 1; nid < cfg::ARRAY_SIZE; ++nid) {
        const auto start = static_cast<size_t>(offsets_key[nid]);
        const auto end = static_cast<size_t>(offsets_key[nid + 1]);
        const size_t slice_size = end - start;

        // An absent node or a node without a path occurrence
        // keeps NODE_UNSEEN_32.
        if (slice_size == 0) continue;

        auto slice_begin = offsets_table.begin() + start;
        const auto slice_end = offsets_table.begin() + end;
        auto right_median = slice_begin + slice_size / 2;

        /*
         * Select the same value as:  sorted_offsets[slice_size / 2]
         * This is the right median for an even number of observations.
         */
        std::nth_element(slice_begin, right_median, slice_end);
        medians[nid] = *right_median;
    }

    return medians;
}


std::unordered_map<uint64_t, uint32_t> compute_edge_weights(const gbwtgraph::GBWTGraph &graph) {
    std::unordered_map<uint64_t, uint32_t> weights;
    weights.reserve(graph.get_edge_count());

    graph.for_each_path_handle([&](const handlegraph::path_handle_t &path) -> bool {
        handlegraph::nid_t previous_nid = 0;
        bool has_previous_node = false;

        graph.for_each_step_in_path(path, [&](const handlegraph::step_handle_t &step) -> bool {
            const handlegraph::handle_t handle = graph.get_handle_of_step(step);
            const handlegraph::nid_t nid = graph.get_id(handle);

            if (has_previous_node) {
                const uint64_t edge_key = detail::pack_node_pair(previous_nid, nid);
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
