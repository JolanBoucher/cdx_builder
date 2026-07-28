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

namespace cdx
{
    void CdxWriter::write_cdx_stream(
        std::ostream& out,
        const std::vector<std::string>& compo_names,
        const std::vector<std::uint32_t>& nid2compo,
        const std::vector<std::uint32_t>& local_idx,
        const std::vector<std::uint32_t>& seq_len,
        const std::uint64_t node_id_offset,
        const std::size_t buffer_size)
    {
        // 1. Validate input vector sizes against expected configuration limits
        if (local_idx.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "local_idx vector doesn't have the expected size\n"
                "local_idx: " + std::to_string(local_idx.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (seq_len.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "sequence length vector doesn't have the expected size\n"
                "sequence length vector: " + std::to_string(seq_len.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (nid2compo.size() != cfg::ARRAY_SIZE) {
            throw std::invalid_argument(
                "the component binder vector doesn't have the expected size\n"
                "component binder vector: " + std::to_string(nid2compo.size()) +
                ", expected: " + std::to_string(cfg::ARRAY_SIZE));
        }
        if (compo_names.size() != cfg::N_COMPO) {
            throw std::invalid_argument(
                "compo_names vector doesn't have the expected size\n"
                "compo_names vector: " + std::to_string(compo_names.size()) +
                ", expected: " + std::to_string(cfg::N_COMPO));
        }

#ifndef NDEBUG
        for (std::size_t i = 0; i < cfg::ARRAY_SIZE; ++i)
        {
            const bool unseen_compo = nid2compo[i] == cfg::NODE_UNSEEN_32;
            const bool unseen_idx = local_idx[i] == cfg::NODE_UNSEEN_32;
            const bool unseen_len = seq_len[i] == cfg::NODE_UNSEEN_32;
            assert(unseen_compo == unseen_idx && unseen_idx == unseen_len);
        }
#endif

        // Pass 1: Compute per-component record counts & 64-bit node ID bounds
        std::vector<ComponentStats> stats(cfg::N_COMPO);
        for (std::size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
            const std::uint32_t cid = nid2compo[i];

            // Explicitly filter out unseen sentinel nodes
            if (cid == cfg::NODE_UNSEEN_32) continue;

            // Validate that component index falls strictly within bounds
            if (cid >= cfg::N_COMPO) {
                throw std::runtime_error(
                    "Invalid component identifier encountered: " + std::to_string(cid));
            }

            // 1-based node ID computation for GBWT convention
            const std::uint64_t nid = static_cast<std::uint64_t>(i) + node_id_offset + 1;
            auto& comp = stats[cid];

            ++comp.n_records;
            comp.nid_min = std::min(comp.nid_min, nid);
            comp.nid_max = std::max(comp.nid_max, nid);
        }

        // Compute starting offsets for each component (CSR / Counting Sort layout)
        std::vector<std::size_t> compo_offsets(cfg::N_COMPO + 1, 0);
        for (std::size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
            compo_offsets[cid + 1] = compo_offsets[cid] + stats[cid].n_records;
        }

        const std::size_t total_records = compo_offsets.back();

        // Invariant: compo_names and nid2compo must already use the
        // final biological component ordering (e.g., chr1, chr2, ..., chrX, chrY, chrM).

        // Pass 2: Build a contiguous 4-byte index mapping array (node_order)
        std::vector<std::uint32_t> node_order(total_records);
        std::vector<std::size_t> current_offsets = compo_offsets;

        for (std::size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
            const std::uint32_t cid = nid2compo[i];
            if (cid == cfg::NODE_UNSEEN_32) continue;

            node_order[current_offsets[cid]++] = static_cast<std::uint32_t>(i);
        }

        // Initialize output serialization buffer
        std::vector<char> buffer(std::max(buffer_size, CdxFormat::COMPONENT_HEADER_SIZE + 256));
        std::size_t cursor = 0;

        auto flush_buffer = [&] {
            if (cursor > 0) {
                out.write(buffer.data(), static_cast<std::streamsize>(cursor));
                if (!out) {
                    throw std::runtime_error("I/O error during CDX stream write.");
                }
                cursor = 0;
            }
        };

        // STEP 1: Pack and write the global file header
        CdxFormat::pack_file_header(buffer.data() + cursor, cfg::N_COMPO);
        cursor += CdxFormat::FILE_HEADER_SIZE;

        // STEP 2: Sequentially write each component and its node records
        for (std::size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
            const auto& name = compo_names[cid];
            const auto& compo_stat = stats[cid];
            const auto name_size = static_cast<std::uint32_t>(name.size());

            const std::size_t header_and_name_bytes = CdxFormat::COMPONENT_HEADER_SIZE + name_size;
            if (cursor + header_and_name_bytes > buffer.size()) {
                flush_buffer();
            }

            // Pack component metadata header
            CdxFormat::pack_component_header(
                buffer.data() + cursor,
                compo_stat.n_records,
                compo_stat.nid_min,
                compo_stat.nid_max,
                name_size
            );
            cursor += CdxFormat::COMPONENT_HEADER_SIZE;

            // Write UTF-8 component name immediately after fixed header
            if (name_size > 0) {
                std::memcpy(buffer.data() + cursor, name.data(), name_size);
                cursor += name_size;
            }

            // Pack and write node records belonging to this component directly from source arrays
            const std::size_t start_idx = compo_offsets[cid];
            const std::size_t end_idx = compo_offsets[cid + 1];

            for (std::size_t pos = start_idx; pos < end_idx; ++pos) {
                const std::uint32_t i = node_order[pos];
                // 1-based node ID computation for GBWT convention
                const std::uint64_t nid = static_cast<std::uint64_t>(i) + node_id_offset + 1;

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

        // Flush remaining byte chunk in the serialization buffer
        flush_buffer();
    }

    // wrapper of write_cdx_stream to write a standard binary file in cdx format
    void CdxWriter::write_cdx_file(
        const std::string &output_path,
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint32_t>& nid2compo,
        const std::vector<std::uint32_t>&local_idx,
        const std::vector<std::uint32_t>& seq_len,
        const uint64_t nid_offset,
        const size_t buffer_size)
    {
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open the output file for writing: " + output_path);
        }

        // the writer engine
        write_cdx_stream(out, compo_names, nid2compo, local_idx, seq_len, nid_offset, buffer_size);
    }

    // wrapper of write_cdx_stream to write a compressed binary file in cdx format with zstd
    void CdxWriter::write_cdx_zstd_file(
        const std::string &output_path,
        const std::vector<std::string> &compo_names,
        const std::vector<std::uint32_t> &nid2compo,
        const std::vector<std::uint32_t> &local_idx,
        const std::vector<std::uint32_t> &seq_len,
        int compression_level,
        std::uint64_t node_id_offset,
        std::size_t buffer_size)
    {

        // 1. writing the intermediary binary stream in memory
        std:: ostringstream raw_stream(std::ios::binary);
        write_cdx_stream(
            raw_stream,
            compo_names,
            nid2compo,
            local_idx,
            seq_len,
            node_id_offset,
            buffer_size
        );

        const std::string raw_data = raw_stream.str();

        // 2. Calculating the maximum compression buffer with libzstd
        const std::size_t max_compressed_size = ZSTD_compressBound(raw_data.size());
        std::vector<char> compressed_buffer(max_compressed_size);

        // 3. Compression in one single data block
        const std::size_t compressed_size = ZSTD_compress(
            compressed_buffer.data(),
            compressed_buffer.size(),
            raw_data.data(),
            raw_data.size(),
            compression_level
        );

        if (ZSTD_isError(compressed_size)) {
            throw std::runtime_error("ZSTD compression failled " +
                std::string(ZSTD_getErrorName(compressed_size))
            );
        }

        // 4. Writing on the disk the compressed data block
        std::ofstream out(output_path, std::ios::binary);
        if (!out) {
            throw std::runtime_error("Failed to open the output file for writing: " + output_path);
        }
        out.write(compressed_buffer.data(), static_cast<std::streamsize>(compressed_size));
        if (!out) {
            throw std::runtime_error("I/O error while writing ZSTD compressed file: " + output_path);
        }
    }

    void TsvWriter::write_tsv(
    std::ostream &out,
    const std::vector<uint32_t> &nid2compo,
    const std::vector<uint32_t> &local_idx,
    const std::vector<uint32_t> &start_pos,
    const std::vector<uint32_t> &seq_len,
    const uint64_t node_id_offset,
    const size_t buffer_size)
{
    // 0. Validation of buffer size
    if (buffer_size < 64) {
        throw std::invalid_argument(
            "buffer_size must be at least 64 bytes, given: " + std::to_string(buffer_size));
    }

    // 1. Validation of the entry vectors
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
        const bool unseen_compo = nid2compo[i] == cfg::NODE_UNSEEN_32;
        const bool unseen_idx   = local_idx[i] == cfg::NODE_UNSEEN_32;
        const bool unseen_start = start_pos[i] == cfg::NODE_UNSEEN_32;
        const bool unseen_len   = seq_len[i] == cfg::NODE_UNSEEN_32;

        assert(unseen_compo == unseen_idx &&
               unseen_idx == unseen_start &&
               unseen_start == unseen_len);
    }
#endif

    // 2. First pass: count elements per component
    std::vector<size_t> compo_counts(cfg::N_COMPO, 0);
    for (size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
        const uint32_t cid = nid2compo[i];
        if (cid == cfg::NODE_UNSEEN_32) continue;
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

    // 3. Second pass: Group nodes by component
    std::vector<uint32_t> node_order(compo_offsets.back());
    std::vector<size_t> current_offsets = compo_offsets;

    for (size_t i = 0; i < cfg::ARRAY_SIZE; ++i) {
        const uint32_t cid = nid2compo[i];
        if (cid == cfg::NODE_UNSEEN_32) continue;
        node_order[current_offsets[cid]++] = static_cast<uint32_t>(i);
    }

    // 4. Third pass: sort records inside each component by local_idx
    for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
        auto begin_iter = node_order.begin() + compo_offsets[cid];
        auto end_iter   = node_order.begin() + compo_offsets[cid + 1];

        std::sort(begin_iter, end_iter, [&local_idx](uint32_t a, uint32_t b) {
            return local_idx[a] < local_idx[b];
        });
    }

    // 5. Write TSV header & check stream state
    out << "component_id\tidx\tnode_id\tstart_pos\tlength\n";
    if (!out) {
        throw std::runtime_error("Failed to write TSV header.");
    }

    // 6. Fast zero-allocation serialization buffer
    std::vector<char> buffer(buffer_size);
    size_t cursor = 0;

    auto flush_buffer = [&] {
        if (cursor > 0) {
            out.write(buffer.data(), static_cast<std::streamsize>(cursor));
            if (!out) {
                throw std::runtime_error("I/O error during TSV stream write.");
            }
            cursor = 0;
        }
    };

    auto append_uint64 = [&](uint64_t value, char trailing_char) {
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

    // 7. Streaming records
    for (size_t cid = 0; cid < cfg::N_COMPO; ++cid) {
        const size_t start_idx = compo_offsets[cid];
        const size_t end_idx   = compo_offsets[cid + 1];

        for (size_t pos = start_idx; pos < end_idx; ++pos) {
            const uint32_t i = node_order[pos];
            const uint64_t nid = static_cast<uint64_t>(i) + node_id_offset;

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
