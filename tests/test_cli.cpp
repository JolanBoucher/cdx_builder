/**
 * @file test_cli.cpp
 * @brief Unit test suite for the CLI argument parsing module (src/cli.hpp/.cpp).
 *
 * Covers, in order:
 *   - prepare_output_filepath (output path derivation / extension replacement)
 *   - parse_args              (CLI11-based argument parsing into a CliArgs struct)
 *
 * NOTE ON TESTABILITY: `parse_args()` used to call `std::exit()` directly on any parsing
 * error, which made it untestable without forking a subprocess (GoogleTest death tests).
 * It has been changed (see cli.cpp) to only `std::exit()` for zero-exit-code CLI11 requests
 * (--help/--version, which are not really errors), and to let genuine parsing errors
 * propagate as `CLI::ParseError` (which derives from std::exception) instead. This lets the
 * error-path tests below use plain `EXPECT_THROW` rather than `EXPECT_EXIT`, and matches
 * main.cpp's existing generic `catch (const std::exception&)` handler, which already
 * catches and reports these safely in production.
 */

#include <gtest/gtest.h>
#include <CLI/CLI.hpp>
#include <string>
#include <vector>
#include "cli.hpp"

namespace {

/**
 * @brief Owns the string storage backing an argc/argv pair, so parse_args() can be called
 *        with a mutable char** without lifetime or const-correctness footguns.
 */
class ArgvHolder {
public:
    explicit ArgvHolder(std::vector<std::string> args) : storage_(std::move(args)) {
        argv_.reserve(storage_.size());
        for (auto& s : storage_) argv_.push_back(s.data());
    }

    [[nodiscard]] int argc() const { return static_cast<int>(argv_.size()); }
    char** argv() { return argv_.data(); }

private:
    std::vector<std::string> storage_;
    std::vector<char*> argv_;
};

} // namespace

// ============================================================================
// 1. DERIVATION DU CHEMIN DE SORTIE (prepare_output_filepath)
// ============================================================================
//
// Invariants / contract:
//   - If `output_filepath` is non-empty, it always takes priority over `input_filepath`,
//     including its own (possibly different) extension.
//   - Only the last '.' in the filename is treated as the extension separator, so multi-dot
//     names like "sample.chr1.gbz" correctly keep "sample.chr1" and only replace ".gbz".
//   - `ext` gets a leading '.' prepended if missing; an empty `ext` removes any extension
//     entirely instead of adding one.
//   - Directory components (parent path) are preserved.
//   - Both paths empty returns "".

TEST(PrepareOutputFilepathTest, DefaultsToInputStemWithCdxExtension) {
    EXPECT_EQ(prepare_output_filepath("", "graph.gbz"), "graph.cdx");
}

TEST(PrepareOutputFilepathTest, InputWithNoExtensionStillGetsOneAppended) {
    EXPECT_EQ(prepare_output_filepath("", "graph"), "graph.cdx");
}

TEST(PrepareOutputFilepathTest, ExplicitOutputTakesPriorityOverInput) {
    EXPECT_EQ(prepare_output_filepath("custom_name.bin", "graph.gbz"), "custom_name.cdx");
}

TEST(PrepareOutputFilepathTest, MultiDotFilenameOnlyStripsTheLastExtension) {
    // Regression guard for a real bug: stripping at the *first* '.' would turn this into
    // "sample.cdx" instead of "sample.chr1.cdx".
    EXPECT_EQ(prepare_output_filepath("", "sample.chr1.gbz"), "sample.chr1.cdx");
}

TEST(PrepareOutputFilepathTest, ParentDirectoryIsPreserved) {
    EXPECT_EQ(prepare_output_filepath("", "/data/graphs/sample.gbz"), "/data/graphs/sample.cdx");
}

TEST(PrepareOutputFilepathTest, ExtensionWithoutLeadingDotGetsOnePrepended) {
    EXPECT_EQ(prepare_output_filepath("", "graph.gbz", "tsv"), "graph.tsv");
}

TEST(PrepareOutputFilepathTest, EmptyExtensionRemovesExtensionEntirely) {
    EXPECT_EQ(prepare_output_filepath("", "graph.gbz", ""), "graph");
}

TEST(PrepareOutputFilepathTest, BothPathsEmptyReturnsEmptyString) {
    EXPECT_EQ(prepare_output_filepath("", ""), "");
}

TEST(PrepareOutputFilepathTest, CustomExtensionIsHonored) {
    EXPECT_EQ(prepare_output_filepath("", "graph.gbz", ".zst"), "graph.zst");
}

// ============================================================================
// 2. ANALYSE DES ARGUMENTS EN LIGNE DE COMMANDE (parse_args)
// ============================================================================
//
// Invariants / contract:
//   - Only `input` is mandatory; every other option has a documented default
//     (max_iterations=100, threshold=0.01, lambda_anchor=0.7, debug=false,
//     compression_level unset).
//   - `-c/--compress` accepts 0 or 1 arguments (`type_size(0, 1)`): given bare, it defaults
//     to level 3 (`default_str("3")`); given with a value, that value is used.
//   - `-c/--compress` and `-d/--debug` are mutually exclusive (`excludes()`).
//   - `-t/--threshold` and `-i/--iteration` must be positive (`CLI::PositiveNumber`).
//   - `-l/--lambda-anchor` must be in `[0.0, 1.0]` (`CLI::Range`).
//   - `-c/--compress`'s value must be in `[1, 22]` (`CLI::Range`).
//   - A missing mandatory `input` argument, or any of the above constraint violations, throws
//     `CLI::ParseError` (see the file-level docstring for why this is a throw and not a
//     process exit).

TEST(ParseArgsTest, OnlyRequiredArgumentUsesDocumentedDefaults) {
    ArgvHolder argv({"cdx_builder", "graph.gbz"});
    const CliArgs args = parse_args(argv.argc(), argv.argv());

    EXPECT_EQ(args.input_file, "graph.gbz");
    EXPECT_EQ(args.max_iterations, 100);
    EXPECT_FLOAT_EQ(args.threshold, 0.01f);
    EXPECT_FLOAT_EQ(args.lambda_anchor, 0.7f);
    EXPECT_EQ(args.output_file, "");
    EXPECT_FALSE(args.compression_level.has_value());
    EXPECT_FALSE(args.debug);
}

TEST(ParseArgsTest, AllLinearizationAndOutputOptionsAreParsed) {
    ArgvHolder argv({
        "cdx_builder", "graph.gbz",
        "-i", "42",
        "-t", "0.005",
        "-l", "0.3",
        "-o", "out.cdx",
    });
    const CliArgs args = parse_args(argv.argc(), argv.argv());

    EXPECT_EQ(args.input_file, "graph.gbz");
    EXPECT_EQ(args.max_iterations, 42);
    EXPECT_FLOAT_EQ(args.threshold, 0.005f);
    EXPECT_FLOAT_EQ(args.lambda_anchor, 0.3f);
    EXPECT_EQ(args.output_file, "out.cdx");
}

TEST(ParseArgsTest, BareCompressFlagDefaultsToLevelThree) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "--compress"});
    const CliArgs args = parse_args(argv.argc(), argv.argv());

    ASSERT_TRUE(args.compression_level.has_value());
    EXPECT_EQ(args.compression_level.value(), 3);
}

TEST(ParseArgsTest, CompressFlagWithExplicitLevelIsHonored) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "--compress", "19"});
    const CliArgs args = parse_args(argv.argc(), argv.argv());

    ASSERT_TRUE(args.compression_level.has_value());
    EXPECT_EQ(args.compression_level.value(), 19);
}

TEST(ParseArgsTest, DebugFlagSetsDebugTrue) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "--debug"});
    const CliArgs args = parse_args(argv.argc(), argv.argv());

    EXPECT_TRUE(args.debug);
}

TEST(ParseArgsTest, MissingRequiredInputThrows) {
    ArgvHolder argv({"cdx_builder"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, NonPositiveThresholdThrows) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "-t", "0"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, NonPositiveIterationCountThrows) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "-i", "-5"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, LambdaAnchorAboveOneThrows) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "-l", "1.5"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, LambdaAnchorBelowZeroThrows) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "-l", "-0.1"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, LambdaAnchorBoundaryValuesAreAccepted) {
    ArgvHolder argv_zero({"cdx_builder", "graph.gbz", "-l", "0.0"});
    EXPECT_NO_THROW(parse_args(argv_zero.argc(), argv_zero.argv()));

    ArgvHolder argv_one({"cdx_builder", "graph.gbz", "-l", "1.0"});
    EXPECT_NO_THROW(parse_args(argv_one.argc(), argv_one.argv()));
}

TEST(ParseArgsTest, CompressionLevelOutOfRangeThrows) {
    ArgvHolder argv_low({"cdx_builder", "graph.gbz", "--compress", "0"});
    EXPECT_THROW(parse_args(argv_low.argc(), argv_low.argv()), CLI::ParseError);

    ArgvHolder argv_high({"cdx_builder", "graph.gbz", "--compress", "23"});
    EXPECT_THROW(parse_args(argv_high.argc(), argv_high.argv()), CLI::ParseError);
}

TEST(ParseArgsTest, DebugAndCompressAreMutuallyExclusive) {
    ArgvHolder argv({"cdx_builder", "graph.gbz", "--debug", "--compress"});
    EXPECT_THROW(parse_args(argv.argc(), argv.argv()), CLI::ParseError);
}
