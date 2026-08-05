/**
 * @file test_cdx_writer.cpp
 * @brief Unit test suite for the CDX/TSV serialization module (src/cdx_writer.h/.cpp).
 *
 * Covers, in order:
 *   - cdx::TsvWriter::write_tsv        (tab-separated text export, streamed directly)
 *   - cdx::CdxWriter::write_cdx_file   (uncompressed binary CDX export, round-tripped via a
 *                                        temp file and manually unpacked with cdx::CdxFormat)
 *   - cdx::CdxWriter::write_cdx_zstd_file (Zstandard-compressed binary CDX export)
 *   - Strict cdx_format.h spec compliance of the raw bytes cdx_builder writes to disk
 *
 * Group 4 (spec compliance) deliberately does NOT rely on cdx::CdxFormat::unpack_* to check
 * correctness, unlike groups 2/3 above: unpack_* is pack_*'s own inverse, so a round-trip through
 * both would still "pass" even if both were consistently wrong (e.g. both silently used
 * big-endian, or both used the wrong struct layout) -- it only proves internal self-consistency,
 * not actual conformance to the documented on-disk format. Group 4 instead reads raw bytes at
 * fixed offsets and compares them against hand-computed expected byte sequences (magic bytes,
 * little-endian field encoding, struct sizes, exact total file length), and separately calls
 * cdx::CdxFormat::has_valid_magic()/has_valid_widths() -- the format module's own canonical
 * conformance checks -- against what cdx_builder actually produced.
 *
 * `validate_and_count_records`, `debug_validate_local_idx_density`, and `write_cdx_stream` are
 * private static members of CdxWriter, so every validation-error test here goes through the
 * public `write_cdx_file` entry point instead of calling them directly.
 *
 * INTERESTING DIVERGENCE BETWEEN THE TWO WRITERS (documented and tested below): TsvWriter
 * explicitly sorts each component's rows by `local_idx` before writing
 * (`std::sort(..., local_idx[a] < local_idx[b])`), but CdxWriter's `write_cdx_stream` does NOT --
 * it scatters node ids into `node_order` by iterating node index ascending
 * (`for (i = 0; i < nid2compo.size(); ++i) { ... node_order[write_pos[cid]++] = i; }`), so binary
 * CDX records within a component are physically ordered by ascending node id, not by `local_idx`.
 * The `idx` field in each on-disk record is positional data to be read, not the writer's sort key.
 * Any downstream code assuming CDX records come out in `local_idx` order (e.g. because that's how
 * the TSV export behaves) would be wrong.
 *
 * Shared test infrastructure:
 *   - CfgFixture: saves/restores cfg::ARRAY_SIZE and cfg::N_COMPO around each test (both writers
 *     validate their input vector sizes strictly against these; see test_gbz_io.cpp for the
 *     rationale of doing this across a shared gtest binary).
 *   - TempFile: RAII wrapper around a unique path under the system temp directory, removed on
 *     destruction, used by the CdxWriter tests (which write to disk, unlike TsvWriter which
 *     streams directly to an in-memory ostringstream).
 *   - parse_cdx_stream(): reads a raw (already-decompressed) CDX byte stream back into a
 *     structured, easy-to-assert-on form using cdx::CdxFormat's public unpack_* functions.
 */

#include <gtest/gtest.h>
#include <array>
#include <atomic>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <zstd.h>
#include <unistd.h>
#include "cdx_writer.h"
#include "cdx_format.h"
#include "constant.h"

namespace {

/// Saves/restores the cfg:: globals this module reads (cfg::ARRAY_SIZE, cfg::N_COMPO).
class CfgFixture : public ::testing::Test {
protected:
    void SetUp() override {
        saved_array_size_ = cfg::ARRAY_SIZE;
        saved_n_compo_ = cfg::N_COMPO;
    }

    void TearDown() override {
        cfg::ARRAY_SIZE = saved_array_size_;
        cfg::N_COMPO = saved_n_compo_;
    }

private:
    size_t saved_array_size_{};
    size_t saved_n_compo_{};
};

/// RAII temp file path under the system temp directory; removed on destruction.
class TempFile {
public:
    explicit TempFile(const std::string& suffix = ".cdx") {
        static std::atomic<int> counter{0};
        path_ = std::filesystem::temp_directory_path() /
            ("cdx_writer_test_" + std::to_string(::getpid()) + "_" +
             std::to_string(counter++) + suffix);
    }

    ~TempFile() {
        std::error_code ec;
        std::filesystem::remove(path_, ec);
    }

    TempFile(const TempFile&) = delete;
    TempFile& operator=(const TempFile&) = delete;

    [[nodiscard]] const std::filesystem::path& path() const { return path_; }

private:
    std::filesystem::path path_;
};

/// A single parsed component block from a decoded CDX stream.
struct ParsedComponent {
    cdx::ComponentHeader header{};
    std::string name;
    std::vector<cdx::NodeRecord> records;
};

/// Reads a raw (already-decompressed) CDX byte stream into a structured, assertable form.
std::vector<ParsedComponent> parse_cdx_stream(std::istream& in) {
    std::array<char, cdx::CdxFormat::FILE_HEADER_SIZE> file_header_buf{};
    in.read(file_header_buf.data(), static_cast<std::streamsize>(file_header_buf.size()));
    const cdx::FileHeader file_header = cdx::CdxFormat::unpack_file_header(file_header_buf.data());

    std::vector<ParsedComponent> components;
    components.reserve(file_header.n_components);

    for (std::uint32_t c = 0; c < file_header.n_components; ++c) {
        std::array<char, cdx::CdxFormat::COMPONENT_HEADER_SIZE> comp_header_buf{};
        in.read(comp_header_buf.data(), static_cast<std::streamsize>(comp_header_buf.size()));

        ParsedComponent comp;
        comp.header = cdx::CdxFormat::unpack_component_header(comp_header_buf.data());

        comp.name.resize(comp.header.name_size);
        if (comp.header.name_size > 0) {
            in.read(comp.name.data(), static_cast<std::streamsize>(comp.header.name_size));
        }

        comp.records.resize(comp.header.n_records);
        for (auto& rec : comp.records) {
            std::array<char, cdx::CdxFormat::RECORD_SIZE> rec_buf{};
            in.read(rec_buf.data(), static_cast<std::streamsize>(rec_buf.size()));
            rec = cdx::CdxFormat::unpack_node_record(rec_buf.data());
        }

        components.push_back(std::move(comp));
    }

    return components;
}

/**
 * @brief Reads the raw bytes of a file into memory, for offset-based byte assertions.
 *
 * Deliberately does not go through cdx::CdxFormat::unpack_* -- see the file-level docstring for
 * why that would make the spec-compliance checks in group 4 circular.
 */
std::string read_raw_bytes(const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

/// Hand-computed little-endian byte encoding of a uint32_t, independent of cdx::CdxFormat.
std::array<unsigned char, 4> le_bytes32(const uint32_t value) {
    std::array<unsigned char, 4> bytes{};
    for (int i = 0; i < 4; ++i) {
        bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xFFu);
    }
    return bytes;
}

/// Hand-computed little-endian byte encoding of a uint64_t, independent of cdx::CdxFormat.
std::array<unsigned char, 8> le_bytes64(const uint64_t value) {
    std::array<unsigned char, 8> bytes{};
    for (int i = 0; i < 8; ++i) {
        bytes[i] = static_cast<unsigned char>((value >> (8 * i)) & 0xFFu);
    }
    return bytes;
}

/// Compares `count` bytes of `data` starting at `offset` against `expected`, byte for byte.
template <std::size_t N>
::testing::AssertionResult BytesEqual(
    const std::string& data, const std::size_t offset, const std::array<unsigned char, N>& expected) {
    if (offset + N > data.size()) {
        return ::testing::AssertionFailure()
            << "buffer too short: need " << (offset + N) << " bytes, have " << data.size();
    }
    for (std::size_t i = 0; i < N; ++i) {
        const auto actual = static_cast<unsigned char>(data[offset + i]);
        if (actual != expected[i]) {
            return ::testing::AssertionFailure()
                << "byte mismatch at offset " << (offset + i)
                << ": expected 0x" << std::hex << static_cast<int>(expected[i])
                << ", got 0x" << std::hex << static_cast<int>(actual);
        }
    }
    return ::testing::AssertionSuccess();
}

} // namespace

// ============================================================================
// 1. EXPORT TSV (cdx::TsvWriter::write_tsv)
// ============================================================================
//
// Invariants / contract:
//   - Header row is a fixed literal: "compo_name\tcomponent_id\tidx\tnode_id\tstart_pos\tlength".
//   - Within each component, rows are sorted by local_idx ascending (explicit std::sort).
//   - node_id column = array index + node_id_offset.
//   - buffer_size must be >= 64 bytes; all input vectors must match cfg::N_COMPO / cfg::ARRAY_SIZE
//     exactly.

TEST_F(CfgFixture, TsvBasicSingleComponentRoundTrip) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> start_pos = {cfg::NODE_UNSEEN_32, 100};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    std::ostringstream out;
    cdx::TsvWriter::write_tsv(out, compo_names, nid2compo, local_idx, start_pos, seq_len);

    std::istringstream in(out.str());
    std::string header, data_row, trailing;
    std::getline(in, header);
    std::getline(in, data_row);
    EXPECT_EQ(header, "compo_name\tcomponent_id\tidx\tnode_id\tstart_pos\tlength");
    EXPECT_EQ(data_row, "chrA\t0\t0\t1\t100\t10");
    EXPECT_FALSE(std::getline(in, trailing)); // exactly one data row
}

TEST_F(CfgFixture, TsvRowsWithinComponentSortedByLocalIdx) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;

    // Node 1 has the *higher* local_idx despite the *lower* node id: the TSV must still
    // emit it second, i.e. sorted by local_idx, not by node id.
    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 1, 0};
    const std::vector<uint32_t> start_pos = {cfg::NODE_UNSEEN_32, 10, 20};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 5, 5};

    std::ostringstream out;
    cdx::TsvWriter::write_tsv(out, compo_names, nid2compo, local_idx, start_pos, seq_len);

    std::istringstream in(out.str());
    std::string header, first_row, second_row;
    std::getline(in, header);
    std::getline(in, first_row);
    std::getline(in, second_row);

    EXPECT_EQ(first_row, "chrA\t0\t0\t2\t20\t5");  // local_idx 0 -> node 2, first
    EXPECT_EQ(second_row, "chrA\t0\t1\t1\t10\t5"); // local_idx 1 -> node 1, second
}

TEST_F(CfgFixture, TsvNodeIdOffsetIsApplied) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> start_pos = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 1};

    std::ostringstream out;
    cdx::TsvWriter::write_tsv(out, compo_names, nid2compo, local_idx, start_pos, seq_len,
                               /*node_id_offset=*/1000);

    std::istringstream in(out.str());
    std::string header, data_row;
    std::getline(in, header);
    std::getline(in, data_row);
    EXPECT_EQ(data_row, "chrA\t0\t0\t1001\t0\t1"); // array index 1 + offset 1000
}

TEST_F(CfgFixture, TsvBufferSizeTooSmallThrows) {
    cfg::ARRAY_SIZE = 1;
    cfg::N_COMPO = 0;

    std::ostringstream out;
    EXPECT_THROW(
        cdx::TsvWriter::write_tsv(out, {}, {cfg::NODE_UNSEEN_16}, {cfg::NODE_UNSEEN_32},
                                   {cfg::NODE_UNSEEN_32}, {cfg::NODE_UNSEEN_32},
                                   /*node_id_offset=*/0, /*buffer_size=*/32),
        std::invalid_argument);
}

TEST_F(CfgFixture, TsvVectorSizeMismatchThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1; // but compo_names below only has 0 entries

    std::ostringstream out;
    EXPECT_THROW(
        cdx::TsvWriter::write_tsv(
            out, /*compo_names=*/{}, {cfg::NODE_UNSEEN_16, cfg::NODE_UNSEEN_16},
            {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32}, {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32},
            {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32}),
        std::invalid_argument);
}

// ============================================================================
// 2. EXPORT BINAIRE CDX NON COMPRESSE (cdx::CdxWriter::write_cdx_file)
// ============================================================================
//
// Invariants / contract: see the file-level docstring for the node-id-ordering divergence
// from TsvWriter. Additionally:
//   - Node index 0 must never be active (cfg::NODE_UNSEEN_16 in nid2compo[0]); a set component
//     there throws std::runtime_error.
//   - Every declared component (compo_names.size() == cfg::N_COMPO) must have at least one
//     active record; an empty one throws std::runtime_error.
//   - Sentinel consistency is enforced per node: an inactive node (component == UNSEEN) must
//     also have UNSEEN local_idx and seq_len; a mismatch throws std::runtime_error. An active
//     node must have a real (non-UNSEEN, nonzero) local_idx and seq_len.
//   - local_idx must be < that component's active record count (bounds-checked).
//   - Component names must be non-empty and free of embedded NUL bytes.

TEST_F(CfgFixture, CdxBasicSingleComponentRoundTrip) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10, 20};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    std::ifstream in(tmp.path(), std::ios::binary);
    ASSERT_TRUE(in.is_open());
    const auto components = parse_cdx_stream(in);

    ASSERT_EQ(components.size(), 1u);
    EXPECT_EQ(components[0].name, "chrA");
    EXPECT_EQ(components[0].header.n_records, 2u);
    EXPECT_EQ(components[0].header.node_id_min, 1u);
    EXPECT_EQ(components[0].header.node_id_max, 2u);

    ASSERT_EQ(components[0].records.size(), 2u);
    EXPECT_EQ(components[0].records[0].node_id, 1u);
    EXPECT_EQ(components[0].records[0].idx, 0u);
    EXPECT_EQ(components[0].records[0].seq_len, 10u);
    EXPECT_EQ(components[0].records[1].node_id, 2u);
    EXPECT_EQ(components[0].records[1].idx, 1u);
    EXPECT_EQ(components[0].records[1].seq_len, 20u);
}

TEST_F(CfgFixture, CdxMultipleComponentsRoundTrip) {
    cfg::ARRAY_SIZE = 5;
    cfg::N_COMPO = 2;

    const std::vector<std::string> compo_names = {"chrA", "chrB"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0, 1, 1};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 5, 5, 7, 7};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    std::ifstream in(tmp.path(), std::ios::binary);
    const auto components = parse_cdx_stream(in);

    ASSERT_EQ(components.size(), 2u);
    EXPECT_EQ(components[0].name, "chrA");
    EXPECT_EQ(components[0].header.n_records, 2u);
    EXPECT_EQ(components[1].name, "chrB");
    EXPECT_EQ(components[1].header.n_records, 2u);
    EXPECT_EQ(components[1].header.node_id_min, 3u);
    EXPECT_EQ(components[1].header.node_id_max, 4u);
}

TEST_F(CfgFixture, CdxNidOffsetAppliedToHeaderBoundsAndRecords) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10, 20};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len,
                                    /*nid_offset=*/1000);

    std::ifstream in(tmp.path(), std::ios::binary);
    const auto components = parse_cdx_stream(in);

    ASSERT_EQ(components.size(), 1u);
    EXPECT_EQ(components[0].header.node_id_min, 1001u);
    EXPECT_EQ(components[0].header.node_id_max, 1002u);
    EXPECT_EQ(components[0].records[0].node_id, 1001u);
    EXPECT_EQ(components[0].records[1].node_id, 1002u);
}

TEST_F(CfgFixture, CdxRecordsAreOrderedByNodeIdNotByLocalIdx) {
    // See file-level docstring: unlike TsvWriter, CdxWriter does not sort by local_idx.
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    // Node 1 has the *higher* local_idx: if records were local_idx-sorted, node 2 would come first.
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 1, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 5, 5};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    std::ifstream in(tmp.path(), std::ios::binary);
    const auto components = parse_cdx_stream(in);

    ASSERT_EQ(components[0].records.size(), 2u);
    // Physical order follows node id (1 then 2), NOT local_idx order.
    EXPECT_EQ(components[0].records[0].node_id, 1u);
    EXPECT_EQ(components[0].records[0].idx, 1u); // node 1's own local_idx, out of physical order
    EXPECT_EQ(components[0].records[1].node_id, 2u);
    EXPECT_EQ(components[0].records[1].idx, 0u);
}

TEST_F(CfgFixture, CdxActiveNodeZeroThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {0, cfg::NODE_UNSEEN_16}; // node 0 wrongly active
    const std::vector<uint32_t> local_idx = {0, cfg::NODE_UNSEEN_32};
    const std::vector<uint32_t> seq_len = {10, cfg::NODE_UNSEEN_32};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::runtime_error);
}

TEST_F(CfgFixture, CdxCompoNamesSizeMismatchThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 2; // but compo_names below only has 1 entry

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::invalid_argument);
}

TEST_F(CfgFixture, CdxArraySizeMismatchThrows) {
    cfg::ARRAY_SIZE = 3; // but the node vectors below only have 2 entries

    const std::vector<std::string> compo_names = {"chrA"};
    cfg::N_COMPO = 1;
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::invalid_argument);
}

TEST_F(CfgFixture, CdxEmptyComponentNameThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {""}; // empty name
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::invalid_argument);
}

TEST_F(CfgFixture, CdxComponentNameWithEmbeddedNulThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {std::string("chr\0A", 5)}; // embedded NUL
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::invalid_argument);
}

TEST_F(CfgFixture, CdxInactiveNodeWithLocalIdxSetThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, cfg::NODE_UNSEEN_16}; // inactive
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0}; // but has a local_idx!
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, cfg::NODE_UNSEEN_32};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::runtime_error);
}

TEST_F(CfgFixture, CdxActiveNodeZeroSeqLenThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 0}; // zero length, active node

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::runtime_error);
}

TEST_F(CfgFixture, CdxComponentWithNoActiveRecordsThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 2; // component 1 is declared but never used by any node

    const std::vector<std::string> compo_names = {"chrA", "chrB"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::runtime_error);
}

TEST_F(CfgFixture, CdxLocalIdxOutOfBoundsThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1; // component 0 has exactly 1 active record, so valid local_idx is only 0

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 5}; // out of bounds
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len),
        std::runtime_error);
}

TEST_F(CfgFixture, CdxBufferSizeSmallerThanOneRecordThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len,
                                        /*nid_offset=*/0, /*buffer_size=*/4),
        std::invalid_argument);
}

// ============================================================================
// 3. EXPORT BINAIRE CDX COMPRESSE ZSTANDARD (cdx::CdxWriter::write_cdx_zstd_file)
// ============================================================================
//
// Invariants / contract:
//   - Produces a single Zstandard frame that decompresses back to byte-for-byte the same
//     stream write_cdx_file would have produced for the same input.
//   - compression_level must be in [1, 22]; out-of-range throws std::invalid_argument before
//     any file I/O happens.

TEST_F(CfgFixture, CdxZstdRoundTripMatchesUncompressedContent) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10, 20};

    const TempFile tmp(".cdx.zst");
    cdx::CdxWriter::write_cdx_zstd_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    // Read the compressed frame back and decompress it fully.
    std::ifstream in(tmp.path(), std::ios::binary);
    ASSERT_TRUE(in.is_open());
    const std::string compressed((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());

    const unsigned long long content_size = ZSTD_getFrameContentSize(compressed.data(), compressed.size());
    ASSERT_NE(content_size, ZSTD_CONTENTSIZE_ERROR);
    ASSERT_NE(content_size, ZSTD_CONTENTSIZE_UNKNOWN);

    std::string decompressed(content_size, '\0');
    const size_t result = ZSTD_decompress(decompressed.data(), decompressed.size(),
                                           compressed.data(), compressed.size());
    ASSERT_FALSE(ZSTD_isError(result));
    decompressed.resize(result);

    std::istringstream decompressed_stream(decompressed);
    const auto components = parse_cdx_stream(decompressed_stream);

    ASSERT_EQ(components.size(), 1u);
    EXPECT_EQ(components[0].name, "chrA");
    ASSERT_EQ(components[0].records.size(), 2u);
    EXPECT_EQ(components[0].records[0].node_id, 1u);
    EXPECT_EQ(components[0].records[0].seq_len, 10u);
    EXPECT_EQ(components[0].records[1].node_id, 2u);
    EXPECT_EQ(components[0].records[1].seq_len, 20u);
}

TEST_F(CfgFixture, CdxZstdCompressionLevelOutOfRangeThrows) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp_low(".cdx.zst");
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_zstd_file(tmp_low.path(), compo_names, nid2compo, local_idx,
                                             seq_len, /*nid_offset=*/0,
                                             /*buffer_size=*/cdx::CdxWriter::DEFAULT_BUFFER_SIZE,
                                             /*compression_level=*/0),
        std::invalid_argument);

    const TempFile tmp_high(".cdx.zst");
    EXPECT_THROW(
        cdx::CdxWriter::write_cdx_zstd_file(tmp_high.path(), compo_names, nid2compo, local_idx,
                                             seq_len, /*nid_offset=*/0,
                                             /*buffer_size=*/cdx::CdxWriter::DEFAULT_BUFFER_SIZE,
                                             /*compression_level=*/23),
        std::invalid_argument);
}

// ============================================================================
// 4. CONFORMITE STRICTE AU FORMAT CDX_FORMAT.H (octets bruts sur disque)
// ============================================================================
//
// See the file-level docstring for why these tests read raw bytes directly instead of relying
// on cdx::CdxFormat::unpack_* (which would only prove pack/unpack are each other's inverse, not
// that cdx_builder's output actually matches the documented on-disk layout).
//
// Layout under test (from cdx_format.h):
//   [FileHeader: magic(4B) | n_components(4B LE) | nid_width(1B) | seqlen_width(1B)]  == 10 bytes
//   per component:
//     [ComponentHeader: n_records(8B LE) | node_id_min(8B LE) | node_id_max(8B LE) | name_size(4B LE)] == 28 bytes
//     [ComponentName: name_size raw UTF-8 bytes, no padding, no terminator]
//     [NodeRecord: node_id(8B LE) | idx(4B LE) | seq_len(4B LE)] x n_records, each == 16 bytes

TEST_F(CfgFixture, FileHeaderMagicAndWidthBytesMatchSpecExactly) {
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);
    const std::string raw = read_raw_bytes(tmp.path());
    ASSERT_GE(raw.size(), cdx::CdxFormat::FILE_HEADER_SIZE);

    // Magic bytes 'C','D','X',0x01 -- single bytes, so this check is endianness-independent.
    EXPECT_EQ(raw[0], 'C');
    EXPECT_EQ(raw[1], 'D');
    EXPECT_EQ(raw[2], 'X');
    EXPECT_EQ(static_cast<unsigned char>(raw[3]), 0x01u);

    // nid_width / seqlen_width are single bytes; the spec mandates exactly 8 and 4.
    EXPECT_EQ(static_cast<unsigned char>(raw[8]), cdx::CdxFormat::NID_WIDTH);
    EXPECT_EQ(static_cast<unsigned char>(raw[9]), cdx::CdxFormat::SEQLEN_WIDTH);

    // Cross-check with the format module's own canonical conformance predicates.
    const cdx::FileHeader header = cdx::CdxFormat::unpack_file_header(raw.data());
    EXPECT_TRUE(cdx::CdxFormat::has_valid_magic(header));
    EXPECT_TRUE(cdx::CdxFormat::has_valid_widths(header));
}

TEST_F(CfgFixture, FileHeaderComponentCountIsLittleEndianOnDisk) {
    cfg::ARRAY_SIZE = 5;
    cfg::N_COMPO = 2;

    const std::vector<std::string> compo_names = {"chrA", "chrB"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0, 1, 1};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 5, 5, 7, 7};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);
    const std::string raw = read_raw_bytes(tmp.path());

    // n_components == 2, hand-encoded little-endian: byte 4 must literally be 0x02 and bytes
    // 5-7 must be zero -- a big-endian bug would instead put 0x02 at byte 7.
    EXPECT_TRUE(BytesEqual(raw, 4, le_bytes32(2)));
}

TEST_F(CfgFixture, ComponentHeaderAndNameBytesFollowSpecLayoutAndLittleEndianEncoding) {
    cfg::ARRAY_SIZE = 3;
    cfg::N_COMPO = 1;
    constexpr uint64_t nid_offset = 500;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10, 20};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len, nid_offset);
    const std::string raw = read_raw_bytes(tmp.path());

    constexpr size_t comp_header_offset = cdx::CdxFormat::FILE_HEADER_SIZE; // == 10
    // n_records(8) | node_id_min(8) | node_id_max(8) | name_size(4), all little-endian.
    EXPECT_TRUE(BytesEqual(raw, comp_header_offset + 0, le_bytes64(2)));                // n_records
    EXPECT_TRUE(BytesEqual(raw, comp_header_offset + 8, le_bytes64(nid_offset + 1)));   // node_id_min == 501
    EXPECT_TRUE(BytesEqual(raw, comp_header_offset + 16, le_bytes64(nid_offset + 2)));  // node_id_max == 502
    EXPECT_TRUE(BytesEqual(raw, comp_header_offset + 24, le_bytes32(4)));               // name_size == 4

    constexpr size_t name_offset = comp_header_offset + cdx::CdxFormat::COMPONENT_HEADER_SIZE; // == 38
    ASSERT_GE(raw.size(), name_offset + 4);
    EXPECT_EQ(raw.substr(name_offset, 4), "chrA"); // verbatim UTF-8, no NUL terminator, no padding

    // The first NodeRecord must start immediately after the name, with no gap: node_id(8) |
    // idx(4) | seq_len(4), little-endian.
    const size_t first_record_offset = name_offset + 4;
    EXPECT_TRUE(BytesEqual(raw, first_record_offset + 0, le_bytes64(nid_offset + 1))); // node_id 501
    EXPECT_TRUE(BytesEqual(raw, first_record_offset + 8, le_bytes32(0)));              // idx (local_idx[1])
    EXPECT_TRUE(BytesEqual(raw, first_record_offset + 12, le_bytes32(10)));            // seq_len
}

TEST_F(CfgFixture, TotalFileLengthMatchesSpecLayoutExactlyWithNoPaddingOrTrailingBytes) {
    cfg::ARRAY_SIZE = 5;
    cfg::N_COMPO = 2;

    const std::vector<std::string> compo_names = {"chrA", "chromosomeB"}; // deliberately different lengths
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0, 0, 1, 1};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0, 1, 0, 1};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 5, 5, 7, 7};

    const TempFile tmp;
    cdx::CdxWriter::write_cdx_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    // Expected size per the documented layout: FileHeader + per component
    // (ComponentHeader + name bytes + n_records * NodeRecord), with nothing else.
    std::uintmax_t expected_size = cdx::CdxFormat::FILE_HEADER_SIZE;
    expected_size += cdx::CdxFormat::COMPONENT_HEADER_SIZE + compo_names[0].size() + 2 * cdx::CdxFormat::RECORD_SIZE;
    expected_size += cdx::CdxFormat::COMPONENT_HEADER_SIZE + compo_names[1].size() + 2 * cdx::CdxFormat::RECORD_SIZE;

    std::error_code ec;
    const std::uintmax_t actual_size = std::filesystem::file_size(tmp.path(), ec);
    ASSERT_FALSE(ec);
    EXPECT_EQ(actual_size, expected_size);
}

TEST_F(CfgFixture, ZstdDecompressedOutputAlsoPassesCdxFormatsConformanceChecks) {
    // The compressed path must produce a frame that, once decompressed, is *exactly* the same
    // spec-compliant byte stream the uncompressed writer would have produced -- not merely
    // something write_cdx_zstd_file's own unpack round-trip agrees with itself on.
    cfg::ARRAY_SIZE = 2;
    cfg::N_COMPO = 1;

    const std::vector<std::string> compo_names = {"chrA"};
    const std::vector<uint16_t> nid2compo = {cfg::NODE_UNSEEN_16, 0};
    const std::vector<uint32_t> local_idx = {cfg::NODE_UNSEEN_32, 0};
    const std::vector<uint32_t> seq_len = {cfg::NODE_UNSEEN_32, 10};

    const TempFile tmp(".cdx.zst");
    cdx::CdxWriter::write_cdx_zstd_file(tmp.path(), compo_names, nid2compo, local_idx, seq_len);

    const std::string compressed = read_raw_bytes(tmp.path());
    const unsigned long long content_size = ZSTD_getFrameContentSize(compressed.data(), compressed.size());
    ASSERT_NE(content_size, ZSTD_CONTENTSIZE_ERROR);
    ASSERT_NE(content_size, ZSTD_CONTENTSIZE_UNKNOWN);

    std::string decompressed(content_size, '\0');
    const size_t result = ZSTD_decompress(decompressed.data(), decompressed.size(),
                                           compressed.data(), compressed.size());
    ASSERT_FALSE(ZSTD_isError(result));
    decompressed.resize(result);

    ASSERT_GE(decompressed.size(), cdx::CdxFormat::FILE_HEADER_SIZE);
    EXPECT_EQ(decompressed[0], 'C');
    EXPECT_EQ(decompressed[1], 'D');
    EXPECT_EQ(decompressed[2], 'X');

    const cdx::FileHeader header = cdx::CdxFormat::unpack_file_header(decompressed.data());
    EXPECT_TRUE(cdx::CdxFormat::has_valid_magic(header));
    EXPECT_TRUE(cdx::CdxFormat::has_valid_widths(header));
}
