/**
 * @file cli.hpp
 * @brief Command-line interface definitions and argument parsing utilities for the CDX builder.
 *
 * This module defines the configuration structure (`CliArgs`) holding pipeline parameters
 * and provides functions to parse command-line options via CLI11, as well as handle output
 * file path generation and extension formatting.
 *
 * @created Jolan on 2026-07-22.
 */

#ifndef CDX_CLI_HPP
#define CDX_CLI_HPP

#include <optional>
#include <string>

/**
 * @brief Container structure for parsed command-line arguments.
 */
struct CliArgs {
    // Input option
    std::string input_file;                 // Path to the input GBZ graph file (positional)

    // Linearization relaxation parameters
    int max_iterations = 100;               // Maximum relaxation iterations (-i, --iteration)
    float threshold = 0.01;                // Convergence threshold (-t, --threshold)
    float lambda_anchor = 0.7;             // Anchor regularization factor [0.0, 1.0] (-l, --lambda-anchor)

    // Output options
    std::string output_file;               // Path to the generated output CDX file (-o, --output)
    std::optional<int> compression_level;  // Zstandard compression level [1-22] if enabled (-c, --compress)
    bool debug = false;                    // Dump TSV to stdout instead of building a CDX index (-d, --debug)
};

/**
 * @brief Parses command-line arguments using CLI11.
 *
 * Validates mandatory inputs, option groups, and mutually exclusive flags,
 * then returns a populated `CliArgs` instance.
 *
 * @param argc Argument count from main.
 * @param argv Argument vector from main.
 * @return CliArgs Populated configuration structure.
 * @throws CLI::ParseError On a genuine parsing error (missing/invalid argument, range
 *         violation, mutually exclusive options, ...). Callers are expected to catch this
 *         (CLI::ParseError derives from std::exception).
 * @note A `--help`/`--version` request (or any other CLI11 request with exit code 0) is
 *       handled internally and terminates the process directly via `std::exit`, matching
 *       standard CLI11 UX -- it does not throw back to the caller.
 */
CliArgs parse_args(int argc, char** argv);

/**
 * @brief Prepares the output file path by ensuring it uses the specified extension.
 *        If `output_filepath` is empty, defaults to using the stem of `input_filepath`.
 *
 * @param output_filepath The user-supplied destination file path (may be empty).
 * @param input_filepath The source input file path used for fallback derivation.
 * @param ext The target file extension (defaults to ".cdx").
 * @return std::string The fully resolved and normalized output file path.
 */
std::string prepare_output_filepath(
    const std::string& output_filepath,
    const std::string& input_filepath,
    std::string ext = ".cdx"
);

#endif