//
// Created by Jolan on 2026-07-27.
//

#include "cdx_writer.h"
#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <zstd.h>
#include "cdx_format.h"
#include "constant.h"
#include <charconv>
#include <string>

namespace cdx {
    // Serializes genomic graph node attributes into an uncompressed binary CDX format stream.
    void CdxWriter::write_cdx_stream(
        std::ostream &out,
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint16_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint32_t> &seq_len,
        const std::uint64_t nid_offset,
        const std::size_t buffer_size) {
        // Validate memory buffer limits against stream size bounds
        if (buffer_size < CdxFormat::RECORD_SIZE) {
            throw std::invalid_argument(
                "CDX output buffer must hold at least one node record (" +
                std::to_string(CdxFormat::RECORD_SIZE) + " bytes)."
            );
        }
        if (buffer_size > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())) {
            throw std::overflow_error("CDX output buffer size exceeds std::streamsize capacity.");
        }

        // Validate vector sizing strictly against configured graph dimensions
        if (compo_names.size() != cfg::N_COMPO) {
            throw std::invalid_argument(
                "compo_names vector size mismatch: " + std::to_string(compo_names.size()) +
                ", expected: " + std::to_string(cfg::N_COMPO));
        }
        if (nid2compo.size() != cfg::ARRAY_SIZE || local_idx.size() != cfg::ARRAY_SIZE || seq_len.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "Input node arrays (nid2compo, local_idx, seq_len) must match global graph array size (" +
                std::to_string(cfg::ARRAY_SIZE) + ").");
        }

        // Ensure node 0 is not active if present in array
        if (!nid2compo.empty() && nid2compo[0] != cfg::NODE_UNSEEN_16) {
            throw std::runtime_error("Active GBZ node ID 0 is invalid.");
        }

        // Lightweight verification & record counting
        const CdxWriteStats stats = validate_and_count_records(
            compo_names,
            nid2compo,
            local_idx,
            seq_len
        );

#ifndef NDEBUG
        // Debug only: Strict permutation & density validation for local_idx
        debug_validate_local_idx_density(
            nid2compo,
            local_idx,
            stats.component_counts
        );
#endif

        const auto &component_counts = stats.component_counts;
        const std::uint64_t total_records_64 = stats.total_records;

        if (total_records_64 > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
            throw std::overflow_error("CDX total record count exceeds system size_t capacity.");
        }
        const auto total_records = static_cast<std::size_t>(total_records_64);

        // CSR Layout: Calculate prefix sums to determine component record boundaries safely
        std::vector<std::size_t> compo_offsets(compo_names.size() + 1, 0);

        for (std::size_t cid = 0; cid < compo_names.size(); ++cid) {
            if (component_counts[cid] > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
                throw std::overflow_error("Component record count exceeds size_t capacity.");
            }

            const auto count = static_cast<std::size_t>(component_counts[cid]);

            if (count > std::numeric_limits<std::size_t>::max() - compo_offsets[cid]) {
                throw std::overflow_error("Component offset overflow detected.");
            }

            compo_offsets[cid + 1] = compo_offsets[cid] + count;
        }

        if (compo_offsets.back() != total_records) {
            throw std::logic_error("Internal inconsistency between component counts and total record count.");
        }

        // Build array mapping CSR positions to original input array indices (node IDs)
        std::vector<std::size_t> node_order(total_records);
        std::vector<std::size_t> write_pos(compo_offsets.begin(), compo_offsets.end() - 1);

        for (std::size_t i = 0; i < nid2compo.size(); ++i) {
            const std::uint16_t cid = nid2compo[i];
            if (cid == cfg::NODE_UNSEEN_16) continue;

            node_order[write_pos[cid]++] = i;
        }

#ifndef NDEBUG
        // Verify component scatter counts in Debug builds without extra allocations
        for (std::size_t cid = 0; cid < compo_names.size(); ++cid) {
            if (write_pos[cid] != compo_offsets[cid + 1]) {
                throw std::logic_error(
                    "[DEBUG] Component scatter count mismatch for CID " + std::to_string(cid)
                );
            }
        }
#endif

        // Dedicated internal buffer reserved strictly for 16-byte fixed node records
        std::vector<char> buffer(buffer_size);
        std::size_t cursor = 0;

        // Helper lambda to flush the internal block buffer to the target output stream
        auto flush_buffer = [&] {
            if (cursor > 0) {
                out.write(buffer.data(), static_cast<std::streamsize>(cursor));
                if (!out) {
                    throw std::runtime_error("I/O error during CDX stream write.");
                }
                cursor = 0;
            }
        };

        // Serialize global file header (fixed size magic & dimensions)
        std::array<char, CdxFormat::FILE_HEADER_SIZE> file_header{};
        CdxFormat::pack_file_header(file_header.data(), static_cast<std::uint32_t>(compo_names.size()));
        out.write(file_header.data(), file_header.size());
        if (!out) {
            throw std::runtime_error("Failed to write CDX global file header.");
        }

        // Component and node serialization loop
        for (std::size_t cid = 0; cid < compo_names.size(); ++cid) {
            const auto &name = compo_names[cid];
            const std::size_t start_idx = compo_offsets[cid];
            const std::size_t end_idx = compo_offsets[cid + 1];

            if (start_idx == end_idx) {
                throw std::logic_error("Empty component reached the CDX serialization stage.");
            }

            const auto n_records = static_cast<std::uint64_t>(end_idx - start_idx);

            // The array index 'i' represents the original GBZ node ID.
            const std::size_t raw_nid_min = node_order[start_idx];
            const std::size_t raw_nid_max = node_order[end_idx - 1];

            if (static_cast<std::uint64_t>(raw_nid_max) > std::numeric_limits<std::uint64_t>::max() - nid_offset) {
                throw std::overflow_error("Node ID overflow after applying nid_offset.");
            }

            const std::uint64_t nid_min = static_cast<std::uint64_t>(raw_nid_min) + nid_offset;
            const std::uint64_t nid_max = static_cast<std::uint64_t>(raw_nid_max) + nid_offset;

            if (name.size() > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())) {
                throw std::overflow_error("Component name length exceeds stream write capacity.");
            }
            const auto name_size = static_cast<std::uint32_t>(name.size());

            // Flush active record buffer prior to writing component header & dynamic string
            flush_buffer();

            // Pack component header using exact canonical layout (28 bytes)
            std::array<char, CdxFormat::COMPONENT_HEADER_SIZE> comp_header{};
            CdxFormat::pack_component_header(
                comp_header.data(),
                n_records,
                nid_min,
                nid_max,
                name_size
            );

            out.write(comp_header.data(), comp_header.size());
            if (!out) {
                throw std::runtime_error("Failed to write CDX component header.");
            }

            // Write variable-length component name string directly following header
            if (name_size > 0) {
                out.write(name.data(), static_cast<std::streamsize>(name.size()));
                if (!out) {
                    throw std::runtime_error("Failed to write CDX component name.");
                }
            }

            // Buffered write loop for 16-byte node records
            for (std::size_t pos = start_idx; pos < end_idx; ++pos) {
                const std::size_t i = node_order[pos];

                if (static_cast<std::uint64_t>(i) > std::numeric_limits<std::uint64_t>::max() - nid_offset) {
                    throw std::overflow_error("Node ID overflow after applying nid_offset.");
                }
                const std::uint64_t nid = static_cast<std::uint64_t>(i) + nid_offset;

                // Flush when the buffer cannot accommodate another 16-byte record
                if (cursor + CdxFormat::RECORD_SIZE > buffer.size()) {
                    flush_buffer();
                }

                CdxFormat::pack_node_record(
                    buffer.data() + cursor,
                    nid,
                    local_idx[i],
                    seq_len[i]
                );
                cursor += CdxFormat::RECORD_SIZE;
            }
        }

        // Flush remaining node records from internal buffer
        flush_buffer();

        // Final stream state validation (Responsibility for out.flush() remains with callers)
        if (!out) {
            throw std::runtime_error("CDX output stream entered a failure state.");
        }
    }

    // Writes an uncompressed binary CDX index file directly to a specified filesystem path.
    void CdxWriter::write_cdx_file(
        const std::filesystem::path &filepath,
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint16_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint32_t> &seq_len,
        const std::uint64_t nid_offset,
        const std::size_t buffer_size) {
        std::ofstream out(filepath, std::ios::binary | std::ios::trunc);

        if (!out.is_open()) {
            throw std::runtime_error("Cannot open CDX output file: " + filepath.string());
        }

        write_cdx_stream(
            out,
            compo_names,
            nid2compo,
            local_idx,
            seq_len,
            nid_offset,
            buffer_size
        );

        out.flush();
        if (!out) {
            throw std::runtime_error("Failed to flush CDX output file: " + filepath.string());
        }

        out.close();
        if (out.fail()) {
            throw std::runtime_error("Failed to close CDX output file: " + filepath.string());
        }
    }

    // Writes a Zstandard-compressed binary CDX index file to disk using full-frame compression.
    void CdxWriter::write_cdx_zstd_file(
        const std::filesystem::path &filepath,
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint16_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint32_t> &seq_len,
        const std::uint64_t nid_offset,
        const std::size_t buffer_size,
        const int compression_level) {
        if (compression_level < 1 || compression_level > 22) {
            throw std::invalid_argument("Zstandard compression level must be between 1 and 22.");
        }

        // 1. Serialize the complete uncompressed CDX stream in memory
        std::ostringstream raw_stream(std::ios::out | std::ios::binary);

        write_cdx_stream(
            raw_stream,
            compo_names,
            nid2compo,
            local_idx,
            seq_len,
            nid_offset,
            buffer_size
        );

        if (!raw_stream) {
            throw std::runtime_error("Failed to generate the uncompressed CDX stream.");
        }
        const std::string raw_data = raw_stream.str();
        if (raw_data.empty()) {
            throw std::runtime_error("Cannot compress an empty CDX stream.");
        }

        // 2. Allocate the maximum required compressed buffer based on Zstandard bounds
        const std::size_t max_compressed_size = ZSTD_compressBound(raw_data.size());
        if (ZSTD_isError(max_compressed_size)) {
            throw std::runtime_error("Failed to calculate the Zstandard compression bound: " +
                                     std::string(ZSTD_getErrorName(max_compressed_size)));
        }

        std::vector<char> compressed_buffer(max_compressed_size);

        // 3. Compress the entire CDX stream as one single Zstandard frame
        const std::size_t compressed_size = ZSTD_compress(
            compressed_buffer.data(),
            compressed_buffer.size(),
            raw_data.data(),
            raw_data.size(),
            compression_level
        );
        if (ZSTD_isError(compressed_size)) {
            throw std::runtime_error(
                "Zstandard compression failed: " + std::string(ZSTD_getErrorName(compressed_size)));
        }

        if (compressed_size > static_cast<std::size_t>(std::numeric_limits<std::streamsize>::max())) {
            throw std::overflow_error("Compressed CDX output exceeds std::streamsize capacity.");
        }

        // 4. Write the single compressed frame to disk
        std::ofstream out(filepath, std::ios::binary | std::ios::trunc);
        if (!out.is_open()) {
            throw std::runtime_error("Failed to open the compressed CDX output file: " + filepath.string());
        }
        out.write(compressed_buffer.data(), static_cast<std::streamsize>(compressed_size));
        if (!out) {
            throw std::runtime_error("I/O error while writing the compressed CDX file: " + filepath.string());
        }
        out.flush();

        if (!out) {
            throw std::runtime_error("Failed to flush the compressed CDX output file: " + filepath.string());
        }
        out.close();

        if (out.fail()) {
            throw std::runtime_error("Failed to close the compressed CDX output file: " + filepath.string());
        }
    }

    // Validates input array dimensions, node sentinel values, and calculates component record counts.
    CdxWriteStats CdxWriter::validate_and_count_records(
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint16_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint32_t> &seq_len) {
        // 1. Check global preconditions and bounds strictly against cfg dimensions
        if (compo_names.empty()) {
            throw std::invalid_argument("Cannot write a CDX file without connected components.");
        }

        if (compo_names.size() != cfg::N_COMPO) {
            throw std::invalid_argument(
                "compo_names vector size mismatch: " + std::to_string(compo_names.size()) +
                ", expected: " + std::to_string(cfg::N_COMPO));
        }

        if (nid2compo.size() != cfg::ARRAY_SIZE || local_idx.size() != cfg::ARRAY_SIZE || seq_len.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "CDX input node vectors must strictly match global ARRAY_SIZE (" + std::to_string(cfg::ARRAY_SIZE) + ").");
        }

        // 2. Validate component names for length limits and prohibited character contents
        for (std::size_t cid = 0; cid < compo_names.size(); ++cid) {
            const std::string &name = compo_names[cid];
            if (name.empty()) {
                throw std::invalid_argument("Connected component " + std::to_string(cid) + " has an empty name.");
            }
            if (name.size() > static_cast<std::size_t>(std::numeric_limits<std::uint32_t>::max())) {
                throw std::overflow_error(
                    "The name of connected component " + std::to_string(cid) + " exceeds uint32_t capacity.");
            }
            if (name.find('\0') != std::string::npos) {
                throw std::invalid_argument("The name of connected component " + std::to_string(cid) +
                                            " contains an embedded NUL byte.");
            }
        }

        CdxWriteStats stats;
        stats.component_counts.assign(compo_names.size(), 0);

        // 3. Pass 1: Verify node sentinel states and aggregate total active records
        for (std::size_t nid = 0; nid < nid2compo.size(); ++nid) {
            const std::uint16_t cid = nid2compo[nid];
            const std::uint32_t idx = local_idx[nid];
            const std::uint32_t length = seq_len[nid];

            // Ensure unseen/inactive nodes maintain consistent sentinel status across all arrays
            if (cid == cfg::NODE_UNSEEN_16) {
                if (idx != cfg::NODE_UNSEEN_32) {
                    throw std::runtime_error(
                        "Inactive node " + std::to_string(nid) + " unexpectedly has a local index.");
                }
                if (length != cfg::NODE_UNSEEN_32) {
                    throw std::runtime_error(
                        "Inactive node " + std::to_string(nid) + " unexpectedly has a sequence length.");
                }
                continue;
            }

            const std::size_t cid_sz = cid;
            if (cid_sz >= compo_names.size()) {
                throw std::runtime_error("Invalid component ID " + std::to_string(cid) +
                                         " for node " + std::to_string(nid) + ".");
            }
            if (idx == cfg::NODE_UNSEEN_32) {
                throw std::runtime_error("Active node " + std::to_string(nid) + " has no local index.");
            }
            if (length == cfg::NODE_UNSEEN_32) {
                throw std::runtime_error("Active node " + std::to_string(nid) + " has no sequence length.");
            }
            if (length == 0) {
                throw std::runtime_error("Active node " + std::to_string(nid) + " has a zero sequence length.");
            }

            if (stats.component_counts[cid_sz] == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Component record count overflow.");
            }
            if (stats.total_records == std::numeric_limits<std::uint64_t>::max()) {
                throw std::overflow_error("Total CDX record count overflow.");
            }

            ++stats.component_counts[cid_sz];
            ++stats.total_records;
        }

        // 4. Ensure every connected component contains at least one record
        for (std::size_t cid = 0; cid < stats.component_counts.size(); ++cid) {
            if (stats.component_counts[cid] == 0) {
                throw std::runtime_error("Connected component " + std::to_string(cid) + " (" + compo_names[cid] + ")"
                                         " contains no records.");
            }
        }

        // 5. Pass 2: Lightweight bounds check verifying local_idx values against component sizes
        for (std::size_t nid = 0; nid < nid2compo.size(); ++nid) {
            const std::uint16_t cid = nid2compo[nid];
            if (cid == cfg::NODE_UNSEEN_16) continue;

            const std::size_t cid_sz = cid;
            if (static_cast<std::uint64_t>(local_idx[nid]) >= stats.component_counts[cid_sz]) {
                throw std::runtime_error("Local index out of bounds for node " + std::to_string(nid) +
                                         ": idx=" + std::to_string(local_idx[nid]) + ", component size=" +
                                         std::to_string(stats.component_counts[cid_sz]) + ".");
            }
        }

        return stats;
    }

#ifndef NDEBUG
    // Debug helper that verifies local index uniqueness, bounds, and contiguous density using bit vectors.
    void CdxWriter::debug_validate_local_idx_density(
        const std::vector<std::uint16_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint64_t> &component_counts) {
        if (nid2compo.size() != local_idx.size()) {
            throw std::logic_error("[DEBUG] nid2compo and local_idx must have identical sizes.");
        }

        std::vector<std::vector<bool> > seen(component_counts.size());

        // Allocate a compact boolean bitmap per component to track index occurrences
        for (std::size_t cid = 0; cid < component_counts.size(); ++cid) {
            const std::uint64_t count = component_counts[cid];
            if (count == 0) {
                throw std::logic_error("[DEBUG] Component " + std::to_string(cid) + " contains no records.");
            }
            if (count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
                throw std::overflow_error("[DEBUG] Component " + std::to_string(cid) + " exceeds size_t capacity.");
            }
            const auto bitmap_size = static_cast<std::size_t>(count);
            if (bitmap_size > seen[cid].max_size()) {
                throw std::length_error("[DEBUG] Local-index bitmap for component " +
                                        std::to_string(cid) + " exceeds vector capacity.");
            }
            seen[cid].assign(bitmap_size, false);
        }

        // Check for local_idx uniqueness and range bounds across active nodes
        for (std::size_t nid = 0; nid < nid2compo.size(); ++nid) {
            const std::uint16_t cid = nid2compo[nid];

            if (cid == cfg::NODE_UNSEEN_16) {
                if (local_idx[nid] != cfg::NODE_UNSEEN_32) {
                    throw std::logic_error("[DEBUG] Inactive node " + std::to_string(nid) +
                                           " unexpectedly has local index " + std::to_string(local_idx[nid]) + ".");
                }
                continue;
            }

            if (static_cast<std::size_t>(cid) >= component_counts.size()) {
                throw std::logic_error("[DEBUG] Invalid component ID " + std::to_string(cid) +
                                       " for node " + std::to_string(nid) + ".");
            }
            const std::uint32_t idx = local_idx[nid];

            if (idx == cfg::NODE_UNSEEN_32) {
                throw std::logic_error("[DEBUG] Active node " + std::to_string(nid) + " has no local index.");
            }
            if (static_cast<std::uint64_t>(idx) >= component_counts[cid]) {
                throw std::logic_error(
                    "[DEBUG] Local index out of bounds for node " + std::to_string(nid) + ": component=" +
                    std::to_string(cid) + ", idx=" + std::to_string(idx) + ", component_size=" +
                    std::to_string(component_counts[cid]) + ".");
            }
            if (seen[cid][idx]) {
                throw std::logic_error("[DEBUG] Duplicate local index in component " + std::to_string(cid) + ": idx=" +
                                       std::to_string(idx) + ", duplicate found at node " + std::to_string(nid) + ".");
            }
            seen[cid][idx] = true;
        }

        // Verify complete contiguous coverage [0, component_size) with no missing indices
        for (std::size_t cid = 0; cid < seen.size(); ++cid) {
            const auto &component_bits = seen[cid];
            for (std::size_t idx = 0; idx < component_bits.size(); ++idx) {
                if (!component_bits[idx]) {
                    throw std::logic_error("[DEBUG] Local-index density violation in component " +
                                           std::to_string(cid) + ": missing idx=" + std::to_string(idx) + ".");
                }
            }
        }
    }

#endif

    // Exports graph component node allocations and sequence offsets to a tab-separated text format (TSV).
    void TsvWriter::write_tsv(
        std::ostream &out,
        const std::vector<std::string> &compo_names,
        const std::vector<uint16_t> &nid2compo,
        const std::vector<uint32_t> &local_idx,
        const std::vector<uint32_t> &start_pos,
        const std::vector<uint32_t> &seq_len,
        const uint64_t node_id_offset,
        const size_t buffer_size) {
        // 0. Ensure memory buffer size can fit minimum field formatting bounds
        if (buffer_size < 64) {
            throw std::invalid_argument(
                "buffer_size must be at least 64 bytes, given: " + std::to_string(buffer_size));
        }

        // 1. Validate matching container lengths against configured graph constants
        if (compo_names.size() != cfg::N_COMPO) {
            throw std::invalid_argument(
                "compo_names vector size mismatch: " + std::to_string(compo_names.size()) +
                ", expected: " + std::to_string(cfg::N_COMPO));
        }
        if (local_idx.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "local_idx vector size mismatch: " + std::to_string(local_idx.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (seq_len.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "seq_len vector size mismatch: " + std::to_string(seq_len.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (nid2compo.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "nid2compo vector size mismatch: " + std::to_string(nid2compo.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (start_pos.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "start_pos vector size mismatch: " + std::to_string(start_pos.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }

#ifndef NDEBUG
        // Sanity check: ensure sentinel status is synchronized across all vectors
        for (size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
            const bool unseen_compo = nid2compo[i] == cfg::NODE_UNSEEN_16;
            const bool unseen_idx = local_idx[i] == cfg::NODE_UNSEEN_32;
            const bool unseen_start = start_pos[i] == cfg::NODE_UNSEEN_32;
            const bool unseen_len = seq_len[i] == cfg::NODE_UNSEEN_32;

            assert(unseen_compo == unseen_idx &&
                unseen_idx == unseen_start &&
                unseen_start == unseen_len);
        }
#endif

        // 2. First pass: count active elements contained in each component
        std::vector<size_t> compo_counts(cfg::N_COMPO, 0);
        for (size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
            const uint16_t cid = nid2compo[i];
            if (cid == cfg::NODE_UNSEEN_16) continue;
            if (cid >= cfg::N_COMPO) {
                throw std::runtime_error("Invalid component identifier " + std::to_string(cid));
            }
            ++compo_counts[cid];
        }

        // CSR prefix sum for offsets
        std::vector<size_t> compo_offsets(cfg::N_COMPO + 1, 0);
        for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
            compo_offsets[cid + 1] = compo_offsets[cid] + compo_counts[cid];
        }

        // 3. Second pass: group node identifiers by component assignment
        std::vector<uint32_t> node_order(compo_offsets.back());
        std::vector<size_t> current_offsets = compo_offsets;

        for (size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
            const uint16_t cid = nid2compo[i];
            if (cid == cfg::NODE_UNSEEN_16) continue;
            node_order[current_offsets[cid]++] = static_cast<uint32_t>(i);
        }

        // 4. Third pass: sort records inside each component by local_idx
        for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
            const auto begin_iter = node_order.begin() + compo_offsets[cid];
            const auto end_iter = node_order.begin() + compo_offsets[cid + 1];

            std::sort(begin_iter, end_iter, [&local_idx](const uint32_t a, const uint32_t b) {
                return local_idx[a] < local_idx[b];
            });
        }

        // 5. Write TSV header & check stream state
        out << "compo_name\tcomponent_id\tidx\tnode_id\tstart_pos\tlength\n";
        if (!out) {
            throw std::runtime_error("Failed to write TSV header.");
        }

        // 6. Fast zero-allocation serialization buffer and helper lambdas
        std::vector<char> buffer(buffer_size);
        size_t cursor = 0;

        // Flushes memory buffer contents directly to stream output
        auto flush_buffer = [&] {
            if (cursor > 0) {
                out.write(buffer.data(), static_cast<std::streamsize>(cursor));
                if (!out) {
                    throw std::runtime_error("I/O error during TSV stream write.");
                }
                cursor = 0;
            }
        };

        // Formats unsigned integer directly into internal buffer without allocations
        auto append_uint64 = [&](const uint64_t value, const char trailing_char) {
            if (cursor + 32 > buffer.size()) {
                flush_buffer();
            }
            auto res = std::to_chars(buffer.data() + cursor, buffer.data() + buffer.size(), value);
            if (res.ec != std::errc{}) {
                flush_buffer();
                res = std::to_chars(buffer.data() + cursor, buffer.data() + buffer.size(), value);
                if (res.ec != std::errc{}) {
                    throw std::runtime_error("std::to_chars failed to serialize value.");
                }
            }
            cursor = static_cast<size_t>(res.ptr - buffer.data());
            buffer[cursor++] = trailing_char;
        };

        // Appends string view content to buffer, falling back to direct stream write if oversized
        auto append_string = [&](const std::string_view str, const char trailing_char) {
            if (cursor + str.size() + 1 > buffer.size()) {
                flush_buffer();
                // If component name size exceeds full buffer, write directly to stream
                if (str.size() + 1 > buffer.size()) {
                    out.write(str.data(), static_cast<std::streamsize>(str.size()));
                    out.put(trailing_char);
                    if (!out) {
                        throw std::runtime_error("I/O error during TSV stream write.");
                    }
                    return;
                }
            }
            std::memcpy(buffer.data() + cursor, str.data(), str.size());
            cursor += str.size();
            buffer[cursor++] = trailing_char;
        };

        // 7. Streaming tab-separated records to output stream
        for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
            const size_t start_idx = compo_offsets[cid];
            const size_t end_idx = compo_offsets[cid + 1];
            const std::string_view name = compo_names[cid];

            for (size_t pos = start_idx; pos < end_idx; ++pos) {
                const uint32_t i = node_order[pos];
                const uint64_t nid = static_cast<uint64_t>(i) + node_id_offset;

                append_string(name, '\t');
                append_uint64(cid, '\t');
                append_uint64(local_idx[i], '\t');
                append_uint64(nid, '\t');
                append_uint64(start_pos[i], '\t');
                append_uint64(seq_len[i], '\n');
            }
        }

        flush_buffer();
    }
}