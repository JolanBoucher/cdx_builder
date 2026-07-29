#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/xpressive/regex_primitives.hpp>
#include <CLI/App.hpp>
#include <gbwtgraph/gbz.h>
#include "constant.h"
#include "gbz_IO.h"
#include "graph_indexing.h"
#include "graph_linearization.h"
#include "cdx_writer.h"
#include "cli.hpp"

//TODO add validation with the validation function to the cdx_writer functions
// /*
// *void validate_cdx(const CdxData& cdx)
// {
//     check_component_order(cdx);
//     check_local_indices(cdx);
//     check_offsets(cdx);
//     check_component_names(cdx);
//     check_coordinate_monotonicity(cdx);
// }
//  */

int main(int argc, char** argv)
{
    try {
        CliArgs args = parse_args(argc, argv);

        auto total_start = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::high_resolution_clock::now();

        std::string gbz_path = args.input_file;
        std::ifstream in(gbz_path, std::ios::binary);
        if (!in) {
            std::cerr << "Erreur : Impossible d'ouvrir " << gbz_path << std::endl;
            return 1;
        }

        gbwtgraph::GBZ gbz;
        gbz.simple_sds_load(in);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cerr << "[INFO] Chargement du fichier GBZ termine en " << elapsed.count() << " s\n" << std::endl;

        //------------------------------------------------------------------------------------

        start = std::chrono::high_resolution_clock::now();

        cfg::NB_NODES = gbz.graph.get_node_count();
        if (cfg::NB_NODES > cfg::NODE_UNSEEN_32) {
            throw std::overflow_error(
                "The number of nodes in this graph exceeds the size allowed of " +
                std::to_string(cfg::NODE_UNSEEN_32)
            );
        }
        const handlegraph::nid_t min_nid = gbz.graph.min_node_id();
        const handlegraph::nid_t max_nid = gbz.graph.max_node_id();

        // calcul de la taille nécessaire des arrays
        cfg::ARRAY_SIZE = static_cast<size_t>(max_nid) + 1;

        // calcul de la densité du graphe
        const auto span = static_cast<double>(max_nid - min_nid + 1);
        const auto density = static_cast<double>(cfg::NB_NODES) / span;

        // Affichage formaté à 3 décimales
        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Found " << cfg::NB_NODES << " nodes (index density: " << density << ") in the graph\n";

        // On peut toujours lire le nombre de nœuds et chemins via gbz.graph
        std::cout << "Paths: " << gbz.graph.get_path_count() << "\n\n";

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Analyse des dimensions terminee en " << elapsed.count() << " s\n" << std::endl;

        //------------------------------------------------------------------------------------
        start = std::chrono::high_resolution_clock::now();

        // load the node length and nid of every node
        // initialize the data_structure at the correct length
        std::vector<uint32_t> n_len(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
        std::vector<uint32_t> parent(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
        std::vector<uint16_t> children(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);

        // retrieve length and nid of every node in the graphe via gbz handle
        gbz.graph.for_each_handle([&](const handlegraph::handle_t& handle) -> bool {
            const handlegraph::nid_t nid = gbz.graph.get_id(handle);
            const size_t node_length = gbz.graph.get_length(handle);

            n_len[nid] = static_cast<uint32_t>(node_length);    // valid node get their length
            parent[nid] = static_cast<uint32_t>(nid);           // valid node start each in their own component
            children[nid] = 0;                                  // valid node start with no child

            return true;
        });

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Initialisation terminee en " << elapsed.count() << " s\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // attribute each node to it's graph's connected component
        std::vector<std::uint16_t> nid2compo = get_graph_components(gbz.graph, parent, children);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Calcul des composantes connexes termine en " << elapsed.count() << " s\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // Associate each component id with it's contig/path name
        size_t max_step_to_check = 0; // 0 for full path validation, or e.g. 100 for partial validation
        std::vector<std::string> component_names = bind_component_names(gbz, nid2compo, max_step_to_check);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Association des noms de composants terminee en " << elapsed.count() << " s\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // evaluate the number of haplotype in this graph
        cfg::N_HAPLO = count_haplotypes(gbz.graph);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Comptage des haplotypes termine en " << elapsed.count() << " s\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // calculate median offset of each node in the graph
        std::vector<uint32_t> offsets_median = compute_nodes_median_offset(gbz.graph);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Calcul des offsets medians termine en " << elapsed.count() << " s\n" << std::endl;

        start = std::chrono::high_resolution_clock::now();

        // calculate weight for each edge in the graph
        std::unordered_map<uint64_t, uint32_t> edges_weight = compute_edge_weights(gbz.graph);

        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cerr << "[INFO] Calcul des poids des aretes termine en " << elapsed.count() << " s\n" << std::endl;

        // return each component name ordered

        //  build CSR table for relaxation
        start = std::chrono::high_resolution_clock::now();
        CSRMatrix matrix = build_csr_matrix(edges_weight);
        end = std::chrono::high_resolution_clock::now();

        elapsed = end - start;
        std::cerr << "[INFO] Construction de la matrice CSR termine en " << elapsed.count() << " s\n" << std::endl;


        //  topology relaxation
        // Paramètres de la relaxation
        float convergence_threshold = args.threshold;
        int max_iterations = args.max_iterations;
        float lambda_anchor = args.lambda_anchor;

        auto start_relax = std::chrono::high_resolution_clock::now();
        auto [relaxed_positions, executed_iterations] = relax_topology(
            matrix,
            n_len,
            offsets_median,
            convergence_threshold,
            max_iterations,
            lambda_anchor
        );
        auto end_relax = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_relax = end_relax - start_relax;

        std::cerr << "[INFO] Relaxation topologique terminee en "
                  << executed_iterations << " iterations ("
                  << elapsed_relax.count() << " s)\n" << std::endl;

        //  midpoint calculation
        auto start_midpoint = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> nodes_midpoint = calculate_midpoint(relaxed_positions, n_len);
        auto end_midpoint = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed_midpoint = end_midpoint - start_midpoint;
        std::cerr << "[INFO] Midpoint calculé en " << elapsed_midpoint.count() << " s\n" << std::endl;

        // building the local index for each component
        auto start_sort = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> idx_table = assign_local_idx(
            nid2compo,
            relaxed_positions,
            nodes_midpoint);
        auto end_sort = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed_sort = end_sort - start_sort;
        std::cerr << "[INFO] Attribution des idx local termine en " << elapsed_sort.count() << " s\n" << std::endl;


        if (args.debug) { // TSV output to stdout
            auto start_tsv_output = std::chrono::high_resolution_clock::now();

            cdx::TsvWriter::write_tsv(
                std::cout,
                component_names,
                nid2compo,
                idx_table,
                relaxed_positions,
                n_len
            );

            auto end_tsv_output = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_tsv = end_tsv_output - start_tsv_output;
            std::cerr << "[INFO] Output written in TSV format in " << elapsed_tsv.count() << " s\n" << std::endl;
        }
        else if (args.compression_level.has_value()) { // Compressed binary CDX output (.cdx.zstd)
            auto start_compressed_output = std::chrono::high_resolution_clock::now();

            const std::string compressed_output_file = prepare_output_filepath(
                args.output_file,
                args.input_file,
                ".cdx.zstd"
            );

            cdx::CdxWriter::write_cdx_zstd_file(
                compressed_output_file,
                component_names,
                nid2compo,
                idx_table,
                n_len,
                args.compression_level.value() // Pass the level
            );

            auto end_compressed_output = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_compression = end_compressed_output - start_compressed_output;
            std::cerr << "[INFO] Output written in compressed CDX format ("
                      << compressed_output_file << ") in "
                      << elapsed_compression.count() << " s\n" << std::endl;
        }
        else { // Standard uncompressed binary CDX output (.cdx)
            auto start_output = std::chrono::high_resolution_clock::now();

            const std::string output_file = prepare_output_filepath(
                args.output_file,
                    args.input_file,
                    ".cdx"
                );

                cdx::CdxWriter::write_cdx_file(
                    output_file,
                    component_names,
                    nid2compo,
                    idx_table,
                    n_len
                );

                auto end_output = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed_output = end_output - start_output;
                std::cerr << "[INFO] Output written in CDX format ("
                          << output_file << ") in "
                          << elapsed_output.count() << " s\n" << std::endl;
            }


        auto total_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_elapsed = total_end - total_start;
        std::cerr << "[INFO] Execution totale terminee en " << total_elapsed.count() << " s" << std::endl;
    }
    catch (const CLI::ParseError& e) {
        // CLI::App exit helper handles printing help message to stdout or stderr
        CLI::App app;
        return app.exit(e);
    }
    catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    return 0;
}
