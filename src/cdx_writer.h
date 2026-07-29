//
// Created by Jolan on 2026-07-27.
//

#ifndef CDX_BUILDER_CDX_WRITER_H
#define CDX_BUILDER_CDX_WRITER_H

#include <string>
#include <vector>

namespace  cdx {
    class CdxWriter  final
    {
    public:
        CdxWriter() = delete;
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024; // 8 MB

        /**
         * @brief Writes a CDX v1 file from array representations of graph nodes.
         *
         * @param output_path Path to the output .cdx file.
         * @param compo_names Vector of UTF-8 component names.
         * @param nid2compo 0-indexed mapping from node id (nid) to component ID.
         * @param local_idx 0-indexed mapping from node id (nid) to local linear position (idx.
         * @param seq_len 0-indexed mapping from node id (nid) to node sequence length.
         * @param nid_offset Offset applied to internal 0-indexed array indices to get 1-based Node IDs (default: 1).
         * @param buffer_size Output write buffer capacity in bytes (default: 8 MB).
         */
        static void write_cdx_file(
            const std::string& output_path,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            std::uint64_t nid_offset = 1ULL,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE
            );

        /**
         * @brief Writing CDX compressed with Zstandard (.cdx.zst).
         */
        static void write_cdx_zstd_file(
            const std::string& output_path,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            int compression_level = 3,
            std::uint64_t node_id_offset = 1ULL,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE
        );


    private:

        /*
         * A writing stream of CDX format v1 that can be passed to any writing format
         */
        static void write_cdx_stream(
            std::ostream& out,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            uint64_t nid_offset = 1ULL,
            size_t buffer_size = DEFAULT_BUFFER_SIZE
            );

        struct ComponentStats
        {
            uint64_t n_records  = 0;
            uint64_t nid_min    = UINT64_MAX;
            uint64_t nid_max    = 0;
        };
    };

    class TsvWriter {
        public:
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024; // 8 MB

        /**
         * @brief Write the cdx to the TSV format (vers std::cout or a dedicated  std::ostream stream).
         * Format : compo_name \t component_id \t idx \t nid \t start_pos \t length \n
         */
        static void write_tsv(
            std::ostream &out,
            const std::vector<std::string>& compo_names,
            const std::vector<uint16_t> &nid2compo,
            const std::vector<uint32_t> &local_idx,
            const std::vector<uint32_t> &start_pos,
            const std::vector<uint32_t> &seq_len,
            uint64_t node_id_offset = 0,
            size_t buffer_size = DEFAULT_BUFFER_SIZE);
    };
};


#endif //CDX_BUILDER_CDX_WRITER_H