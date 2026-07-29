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
#include <iomanip>
#include <stdexcept>
#include <filesystem>

//TODO pass the error message from everywhere to the main and deal with it properly
//TODO give proper commentary and docstring to each function (almost finish)
//TODO check the cMake for the weird warning


// Helper to format timing outputs consistently
inline void print_STEP_time(const std::string& STEP_name, double seconds) {
    std::cerr << "  - " << std::left << std::setw(38) << STEP_name
              << " Completed in " << std::right << std::setw(7) << std::fixed
              << std::setprecision(3) << seconds << " s\n";
}


int main(int argc, char** argv)
{
    try {
        CliArgs args = parse_args(argc, argv);

        const auto total_start = std::chrono::high_resolution_clock::now();

        std::cerr << "======================================================================\n";
        std::cerr << "                         CDX BUILDER PIPELINE                         \n";
        std::cerr << "======================================================================\n";

        //------------------------------------------------------------------------------------
        // 1. Graph Loading & Memory Allocation
        //------------------------------------------------------------------------------------
        std::cerr << "[STEP 1/6] Graph Loading & Memory Allocation\n";
        auto start = std::chrono::high_resolution_clock::now();

        std::string gbz_path = args.input_file;
        std::ifstream in(gbz_path, std::ios::binary);
        if (!in) {
            std::cerr << "[ERROR] Unable to open input GBZ file: " << gbz_path << std::endl;
            return EXIT_FAILURE;
        }

        gbwtgraph::GBZ gbz;
        gbz.simple_sds_load(in);

        cfg::NB_NODES = gbz.graph.get_node_count();
        if (cfg::NB_NODES > cfg::NODE_UNSEEN_32) {
            throw std::overflow_error(
                "Graph node count (" + std::to_string(cfg::NB_NODES) +
                ") exceeds max supported threshold (" + std::to_string(cfg::NODE_UNSEEN_32) + ")"
            );
        }

        const handlegraph::nid_t max_nid = gbz.graph.max_node_id();
        cfg::ARRAY_SIZE = static_cast<size_t>(max_nid) + 1;

        // Allocate buffers based on graph size
        std::vector<uint32_t> n_len(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
        std::vector<uint32_t> parent(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
        std::vector<uint16_t> children(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);

        gbz.graph.for_each_handle([&](const handlegraph::handle_t& handle) -> bool {
            const handlegraph::nid_t nid = gbz.graph.get_id(handle);
            const size_t node_length = gbz.graph.get_length(handle);

            n_len[nid] = static_cast<uint32_t>(node_length);
            parent[nid] = static_cast<uint32_t>(nid);
            children[nid] = 0;

            return true;
        });

        auto end = std::chrono::high_resolution_clock::now();

        const std::string load_summary = "Loaded GBZ graph (" + std::to_string(cfg::NB_NODES) +
                                         " nodes across " + std::to_string(gbz.graph.get_path_count()) + " paths)";
        print_STEP_time(load_summary, std::chrono::duration<double>(end - start).count());

        //------------------------------------------------------------------------------------
        // 2. Connected Components & Path Binding
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 2/6] Connected Components & Path Binding\n";

        start = std::chrono::high_resolution_clock::now();
        std::vector<std::uint16_t> nid2compo = get_graph_components(gbz.graph, parent, children);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("Connected components identified", std::chrono::duration<double>(end - start).count());

        start = std::chrono::high_resolution_clock::now();
        size_t max_step_to_check = 0; // 0 for full path validation
        std::vector<std::string> component_names = bind_component_names(gbz, nid2compo, max_step_to_check);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("Component path names bound", std::chrono::duration<double>(end - start).count());

        start = std::chrono::high_resolution_clock::now();
        cfg::N_HAPLO = count_haplotypes(gbz.graph);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("Haplotypes evaluated", std::chrono::duration<double>(end - start).count());

        //------------------------------------------------------------------------------------
        // 3. Metric Calculation & CSR Construction
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 3/6] Metric Calculation & CSR Construction\n";

        start = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> offsets_median = compute_nodes_median_offset(gbz.graph);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("Node median offsets computed", std::chrono::duration<double>(end - start).count());

        start = std::chrono::high_resolution_clock::now();
        std::unordered_map<uint64_t, uint32_t> edges_weight = compute_edge_weights(gbz.graph);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("Edge weights calculated", std::chrono::duration<double>(end - start).count());

        start = std::chrono::high_resolution_clock::now();
        CSRMatrix matrix = build_csr_matrix(edges_weight);
        end = std::chrono::high_resolution_clock::now();
        print_STEP_time("CSR matrix constructed", std::chrono::duration<double>(end - start).count());

        //------------------------------------------------------------------------------------
        // 4. Topological Relaxation
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 4/6] Topological Relaxation\n";

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
        print_STEP_time("Topology relaxed (" + std::to_string(executed_iterations) + " iterations)",
                         std::chrono::duration<double>(end_relax - start_relax).count());

        //------------------------------------------------------------------------------------
        // 5. Local Indexing & Midpoint Calculation
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 5/6] Local Indexing & Midpoint Calculation\n";

        auto start_midpoint = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> nodes_midpoint = calculate_midpoint(relaxed_positions, n_len);
        auto end_midpoint = std::chrono::high_resolution_clock::now();
        print_STEP_time("Node midpoints calculated", std::chrono::duration<double>(end_midpoint - start_midpoint).count());

        auto start_sort = std::chrono::high_resolution_clock::now();
        std::vector<uint32_t> idx_table = assign_local_idx(
            nid2compo,
            relaxed_positions,
            nodes_midpoint
        );
        auto end_sort = std::chrono::high_resolution_clock::now();
        print_STEP_time("Local component indexing assigned", std::chrono::duration<double>(end_sort - start_sort).count());

        //------------------------------------------------------------------------------------
        // 6. Output Serialization
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 6/6] Output Serialization\n";

        if (args.debug) {
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
            print_STEP_time("TSV output written to stdout", std::chrono::duration<double>(end_tsv_output - start_tsv_output).count());
        }
        else if (args.compression_level.has_value()) {
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
                args.compression_level.value()
            );

            auto end_compressed_output = std::chrono::high_resolution_clock::now();
            print_STEP_time("Compressed CDX output saved (" + compressed_output_file + ")",
                             std::chrono::duration<double>(end_compressed_output - start_compressed_output).count());
        }
        else {
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
            print_STEP_time("Binary CDX output saved (" + output_file + ")",
                             std::chrono::duration<double>(end_output - start_output).count());
        }

        //------------------------------------------------------------------------------------
        // Total Execution Summary
        //------------------------------------------------------------------------------------
        auto total_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_elapsed = total_end - total_start;

        std::cerr << "\n======================================================================\n";
        std::cerr << " [SUCCESS] Pipeline execution finished in "
                  << std::fixed << std::setprecision(3) << total_elapsed.count() << " s\n";
        std::cerr << "======================================================================\n";

    }
    catch (const CLI::ParseError& e) {
        CLI::App app;
        return app.exit(e);
    }
    catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
