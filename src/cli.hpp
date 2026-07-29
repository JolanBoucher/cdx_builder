#ifndef CDX_CLI_HPP
#define CDX_CLI_HPP

#include <optional>
#include <string>

/**
 * @brief Container structure for parsed command-line arguments.
 */
struct CliArgs {
    // Input option
    std::string input_file;    ///< Path to the input GBZ graph file (-x, --input)

    // Linearization relaxation parameters
    int max_iterations = 100;    ///< Maximum relaxation iterations (-i, --iteration)
    double threshold = 0.01;     ///< Convergence threshold (-t, --threshold)
    double lambda_anchor = 0.7;  ///< Anchor regularization factor [0.0, 1.0] (-l, --lambda-anchor)

    // Output options
    std::string output_file;               ///< Path to the generated output CDX file (-o, --output)
    std::optional<int> compression_level;  ///< Zstandard compression level [1-22] if enabled (-c, --compress)
    bool debug = false;                    ///< Dump TSV to stdout instead of building a CDX index (-d, --debug)
};

/**
 * @brief Parses command-line arguments using CLI11.
 * @throws CLI::ParseError On parsing error or help request.
 */
CliArgs parse_args(int argc, char** argv);

/**
 * @brief Prepares the output file path by ensuring it uses the specified extension.
 *        If `output_filepath` is empty, defaults to using the stem of `input_filepath`.
 */
std::string prepare_output_filepath(
    const std::string& output_filepath,
    const std::string& input_filepath,
    std::string ext = ".cdx"
);

#endif // CDX_CLI_HPP