/**
 * @file main.cpp
 * @brief Entry point for the CDX Index Builder executable.
 *
 * This file contains the main program logic for parsing command-line
 * arguments, managing top-level execution flow, catching runtime errors,
 * and invoking the CDX graph processing pipeline (`run_pipeline`).
 *
 * @details
 * The overall execution flow follows three phases:
 * 1. **CLI Parsing:** Parses input parameters (GBZ file path, output destination,
 *    topological relaxation thresholds, compression flags) into a `CliArgs` struct.
 * 2. **Pipeline Execution:** Passes the validated options to `run_pipeline` which
 *    executes the six core index building stages.
 * 3. **Error Handling & Exit:** Intercepts standard exceptions (`std::exception`),
 *    logs user-friendly error messages to `std::cerr`, and returns appropriate
 *    system exit codes.
 *
 * @author Jolan
 * @date 2026
 */

#include "cdx_writer.h"
#include "cli.hpp"
#include "constant.h"
#include "gbz_IO.h"
#include "graph_indexing.h"
#include "graph_linearization.h"

#include <chrono>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cstddef>

#include <gbwtgraph/gbz.h>

/**
 * @brief RAII timer measuring execution duration with steady_clock and exception tracking.
 */
class ScopedTimer {
public:
    using Clock = std::chrono::steady_clock;

    /**
     * @brief Construct timer and capture current clock time and active exception count.
     * @param step_name Label displayed upon timer completion or failure.
     */
    explicit ScopedTimer(std::string step_name)
        : name_(std::move(step_name)),
          start_(Clock::now()),
          uncaught_on_entry_(std::uncaught_exceptions()) {}
    /**
     * @brief Dynamically updates the step label displayed upon timer destruction.
     * @param new_name The new string description for the timed step.
     */
    void update_name(std::string new_name) {
        name_ = std::move(new_name);
    }

    /**
     * @brief Destructor that computes elapsed time and logs completion or failure status.
     * @note Marked noexcept to prevent throwing during unwinding.
     */
    ~ScopedTimer() noexcept {
        try {
            const auto end = Clock::now();
            const double seconds = std::chrono::duration<double>(end - start_).count();
            const bool failed = std::uncaught_exceptions() > uncaught_on_entry_;
            print_step_time(name_, seconds, failed);
        } catch (...) {
            // Guarantee noexcept behavior in destructor
        }
    }

    // Prevent copy/move to avoid duplicate measurement prints
    ScopedTimer(const ScopedTimer&) = delete;
    ScopedTimer& operator=(const ScopedTimer&) = delete;
    ScopedTimer(ScopedTimer&&) = delete;
    ScopedTimer& operator=(ScopedTimer&&) = delete;

private:
    std::string name_;        ///< Name of the pipeline step being measured.
    Clock::time_point start_; ///< Start timestamp recorded on construction.
    int uncaught_on_entry_;   ///< Number of uncaught exceptions present at construction.

    /**
     * @brief Formats and prints the timing output to std::cerr.
     * @param step_name Label of the executed step.
     * @param seconds Elapsed duration in seconds.
     * @param failed True if step failed due to an uncaught exception.
     */
    static void print_step_time(const std::string& step_name, const double seconds, const bool failed) {
        std::cerr << "  - " << std::left << std::setw(50) << step_name
                  << (failed ? " Failed after   " : " Completed in ")
                  << std::right << std::setw(7) << std::fixed
                  << std::setprecision(4) << seconds << " s\n";
    }
};

/**
 * @brief Orchestrates the complete end-to-end execution of the CDX index building pipeline.
 *
 * This function acts as the primary pipeline manager, driving the transformation of a
 * pangenome graph into a serialized CDX index across six distinct steps:
 *
 * 1. **Graph Loading & Allocation:** Loads the input GBZ graph file, validates node counts/IDs,
 *    and initializes sequence length (`n_len`) buffers along with component discovery arrays.
 * 2. **Component & Path Binding:** Computes connected components across the graph topology,
 *    binds path metadata names to component IDs, and counts total graph haplotypes.
 * 3. **Metric Calculation & CSR Construction:** Calculates node median offsets and
 *    path-weighted edge frequencies to construct a high-performance Sparse Row (CSR) matrix.
 * 4. **Topological Relaxation:** Solves node spatial layouts via iterative topological
 *    relaxation using node lengths and CSR connectivity.
 * 5. **Local Indexing:** Computes geometric node midpoints and assigns continuous local
 *    coordinate indexes across graph components.
 * 6. **Output Serialization:** Serializes index tables, component metadata, and node position
 *    mappings to standard output (debug TSV), uncompressed binary `.cdx`, or Zstandard compressed
 *    `.cdx.zst` format.
 *
 * @param args Evaluated command-line configuration arguments specifying input/output file paths,
 *             relaxation convergence criteria (`threshold`, `max_iterations`, `lambda_anchor`),
 *             and output format options (`debug`, `compression_level`).
 *
 * @pre `args.input_file` must point to a readable, valid binary GBZ pangenome graph.
 *
 * @throws std::runtime_error If the GBZ file cannot be opened, the loaded graph is empty
 *                            (0 nodes), or a node ID returns a negative value.
 * @throws std::overflow_error If graph node count exceeds `cfg::NODE_UNSEEN_32`, node ID
 *                             exceeds addressable system memory range, or sequence length
 *                             exceeds maximum supported `uint32_t` capacity.
 *
 * @note **Resource Management:** This function relies strictly on RAII scope isolation
 *       and explicit buffer deallocations (e.g., `gbz = gbwtgraph::GBZ();`, `.clear()`,
 *       and vector swapping) to release intermediate data structures immediately after
 *       their dependent step completes, minimizing peak Resident Set Size (RSS).
 */
static void run_pipeline(const CliArgs& args) {
    using Clock = std::chrono::steady_clock;
    const auto total_start = Clock::now();

    std::cerr << "=============================================================================\n";
    std::cerr << "                                 CDX BUILDER                                 \n";
    std::cerr << "=============================================================================\n";

    //------------------------------------------------------------------------------------
    // 1. Graph Loading & Memory Allocation
    //------------------------------------------------------------------------------------
    std::cerr << "[STEP 1/6] Graph Loading & Memory Allocation\n";

    const std::string gbz_path = args.input_file;
    gbwtgraph::GBZ gbz;

    {
        ScopedTimer t("Loading GBZ graph...");

        std::ifstream in(gbz_path, std::ios::binary);
        if (!in) {
            throw std::runtime_error("Unable to open input GBZ file: " + gbz_path);
        }
        gbz.simple_sds_load(in);

        // Initial validation and metrics
        cfg::NB_NODES = gbz.graph.get_node_count();
        if (cfg::NB_NODES == 0) {
            throw std::runtime_error("The input GBZ graph contains no nodes.");
        }

        const std::size_t m_paths = gbz.index.sequences(); // or gbz.header.paths

        // Update label formatted
        t.update_name("GBZ graph loaded: " + std::to_string(cfg::NB_NODES) +
                      " nodes on " + std::to_string(m_paths) + " paths");
    }

    // Threshold validation
    if (cfg::NB_NODES > cfg::NODE_UNSEEN_32) {
        throw std::overflow_error(
            "Graph node count (" + std::to_string(cfg::NB_NODES) +
            ") exceeds max supported threshold (" + std::to_string(cfg::NODE_UNSEEN_32) + ")"
        );
    }

    const handlegraph::nid_t max_nid = gbz.graph.max_node_id();
    if (max_nid < 0) {
        throw std::runtime_error("GBZ graph returned a negative maximum node ID.");
    }

    const auto max_nid_u64 = static_cast<std::uint64_t>(max_nid);
    if (max_nid_u64 >= static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::overflow_error("Maximum node ID exceeds addressable array capacity.");
    }

    cfg::ARRAY_SIZE = static_cast<std::size_t>(max_nid_u64) + 1;

    std::vector<std::uint32_t> n_len(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
    std::vector<std::uint16_t> nid2compo;

    // Scope isolation for graph setup arrays ('parent' and 'children' deallocated early)
    {
        std::vector<std::uint32_t> parent(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_32);
        std::vector<std::uint16_t> children(cfg::ARRAY_SIZE, cfg::NODE_UNSEEN_16);

        {
            ScopedTimer t("Node metadata arrays initialized");
            gbz.graph.for_each_handle([&](const handlegraph::handle_t& handle) -> bool {
                const handlegraph::nid_t nid = gbz.graph.get_id(handle);
                const std::size_t node_length = gbz.graph.get_length(handle);

                if (node_length >= cfg::NODE_UNSEEN_32) {
                    throw std::overflow_error(
                        "Node " + std::to_string(nid) +
                        " sequence length exceeds supported uint32_t range."
                    );
                }

                n_len[nid] = static_cast<std::uint32_t>(node_length);
                parent[nid] = static_cast<std::uint32_t>(nid);
                children[nid] = 0;

                return true;
            });
        }

        //------------------------------------------------------------------------------------
        // 2. Connected Components & Path Binding
        //------------------------------------------------------------------------------------
        std::cerr << "\n[STEP 2/6] Connected Components & Path Binding\n";

        {
            ScopedTimer t("Connected components identified");
            nid2compo = get_graph_components(gbz.graph, parent, children);

            t.update_name("Identified " + std::to_string(cfg::N_COMPO) + " connected components");
        }
    } // 'parent' and 'children' memory naturally freed here upon scope exit

    std::vector<std::string> component_names;
    {
        ScopedTimer t("Component path names binding");
        constexpr std::size_t max_step_to_check = 0; // 0 for full validation
        component_names = bind_component_names(gbz, nid2compo, max_step_to_check);
        t.update_name("Bound " + std::to_string(component_names.size()) + " components path names");
    }

    {
        ScopedTimer t("Haplotypes evaluation");
        cfg::N_HAPLO = count_haplotypes(gbz.graph);
        t.update_name("Found " + std::to_string(cfg::N_HAPLO) + " haplotypes in the graph");
    }

    //------------------------------------------------------------------------------------
    // 3. Metric Calculation & CSR Construction
    //------------------------------------------------------------------------------------
    std::cerr << "\n[STEP 3/6] Metric Calculation & CSR Construction\n";

    std::vector<std::uint32_t> offsets_median;
    {
        ScopedTimer t("Node median offsets computed");
        offsets_median = compute_nodes_median_offset(gbz.graph);
    }

    std::unordered_map<std::uint64_t, std::uint32_t> edges_weight;
    {
        ScopedTimer t("Edge weights calculated");
        edges_weight = compute_edge_weights(gbz.graph);
    }

    // Deallocate GBZ graph object prior to matrix construction
    gbz = gbwtgraph::GBZ();

    CSRMatrix matrix;
    {
        ScopedTimer t("CSR matrix constructed");
        matrix = build_csr_matrix(edges_weight);
    }

    // Force memory liberation of intermediate hash table buckets
    std::unordered_map<std::uint64_t, std::uint32_t>().swap(edges_weight);

    //------------------------------------------------------------------------------------
    // 4. Topological Relaxation
    //------------------------------------------------------------------------------------
    std::cerr << "\n[STEP 4/6] Topological Relaxation\n";

    const float convergence_threshold = args.threshold;
    const int max_iterations = args.max_iterations;
    const float lambda_anchor = args.lambda_anchor;

    std::vector<std::uint32_t> relaxed_positions;

    {
        ScopedTimer t("Relaxing topology...");

        auto [pos, iters] = relax_topology(
            matrix,
            n_len,
            offsets_median,
            convergence_threshold,
            max_iterations,
            lambda_anchor
        );
        relaxed_positions = std::move(pos);

        t.update_name("Topology relaxed in " + std::to_string(iters) + " iterations");
    }

    // Deallocate relaxation intermediates
    matrix.clear();
    std::vector<std::uint32_t>().swap(offsets_median);

    //------------------------------------------------------------------------------------
    // 5. Local Indexing & Midpoint Calculation
    //------------------------------------------------------------------------------------
    std::cerr << "\n[STEP 5/6] Local Indexing & Midpoint Calculation\n";

    std::vector<std::uint32_t> idx_table;
    {
        std::vector<std::uint32_t> nodes_midpoint;
        {
            ScopedTimer t("Node midpoints calculated");
            nodes_midpoint = calculate_midpoint(relaxed_positions, n_len);
        }

        {
            ScopedTimer t("Local component indexing assigned");
            idx_table = assign_local_idx(
                nid2compo,
                relaxed_positions,
                nodes_midpoint
            );
        }
    } // 'nodes_midpoint' naturally freed here upon scope exit

    //------------------------------------------------------------------------------------
    // 6. Output Serialization
    //------------------------------------------------------------------------------------
    std::cerr << "\n[STEP 6/6] Output Serialization\n";

    if (args.debug) {
        ScopedTimer t("TSV output written to stdout");
        cdx::TsvWriter::write_tsv(
            std::cout,
            component_names,
            nid2compo,
            idx_table,
            relaxed_positions,
            n_len
        );
    }
    else if (args.compression_level.has_value()) {
        const std::string compressed_output_file = prepare_output_filepath(
            args.output_file,
            args.input_file,
            ".cdx.zst"
        );

        ScopedTimer t("Compressed CDX output saved (" + compressed_output_file + ")");
        cdx::CdxWriter::write_cdx_zstd_file(
            compressed_output_file,
            component_names,
            nid2compo,
            idx_table,
            n_len,
            0, // nid_offset
            cdx::CdxWriter::DEFAULT_BUFFER_SIZE,
            args.compression_level.value()
        );
    }
    else {
        const std::string output_file = prepare_output_filepath(
            args.output_file,
            args.input_file,
            ".cdx"
        );

        ScopedTimer t("Binary CDX output saved (" + output_file + ")");
        cdx::CdxWriter::write_cdx_file(
            output_file,
            component_names,
            nid2compo,
            idx_table,
            n_len
        );
    }

    //------------------------------------------------------------------------------------
    // Total Execution Summary
    //------------------------------------------------------------------------------------
    const auto total_end = Clock::now();
    const std::chrono::duration<double> total_elapsed = total_end - total_start;

    std::cerr << "\n=============================================================================\n";
    std::cerr << " [SUCCESS] Execution finished in "
              << std::fixed << std::setprecision(3) << total_elapsed.count() << " s\n";
    std::cerr << "=============================================================================\n";
}

int main(const int argc, char** argv)
{
    try {
        const CliArgs args = parse_args(argc, argv);
        run_pipeline(args);
        return EXIT_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "\n[FATAL ERROR] " << e.what() << '\n';
        return EXIT_FAILURE;
    }
    catch (...) {
        std::cerr << "\n[FATAL ERROR] Unknown non-standard error occurred.\n";
        return EXIT_FAILURE;
    }
}