/**
 * @file cdx_writer.h
 * @brief Binary (CDX) and plain-text (TSV) file serialization utilities for graph genomic nodes.
 *
 * This module provides serialization implementations for writing node-level
 * graph alignments, component mappings, and sequence metrics into formatted output streams.
 * It supports uncompressed binary CDX, Zstandard-compressed binary (.cdx.zst), and tab-separated
 * values (TSV) formats, incorporating debug-time topological density and index permutation validation.
 */

#ifndef CDX_BUILDER_CDX_WRITER_H
#define CDX_BUILDER_CDX_WRITER_H

#include <cstdint>
#include <filesystem>
#include <iosfwd>
#include <string>
#include <vector>
#include <cstdint>

namespace cdx {

    /**
     * @brief Aggregated record metrics generated during the CDX verification pass.
     */
    struct CdxWriteStats
    {
        std::vector<std::uint64_t> component_counts;    // Number of active nodes per component ID.
        std::uint64_t total_records = 0;                // Total count of non-sentinel nodes across all components.
    };

    /**
     * @brief Utility class providing static methods to serialize graph node attributes to CDX format streams and files.
     */
    class CdxWriter final
    {
    public:
        CdxWriter() = delete;
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024; // 8 MB

        /**
         * @brief Serializes graph node attributes into an uncompressed binary CDX index file at the given path.
         *
         * @param filepath Target filesystem path for the output .cdx file.
         * @param compo_names Vector of UTF-8 component names (must match cfg::N_COMPO).
         * @param nid2compo Mapping from node ID to component ID (size must match cfg::ARRAY_SIZE).
         * @param local_idx Mapping from node ID to component-local position (size must match cfg::ARRAY_SIZE).
         * @param seq_len Mapping from node ID to sequence length in base pairs (size must match cfg::ARRAY_SIZE).
         * @param nid_offset Offset added to array indices to compute output Node IDs (default: 0).
         * @param buffer_size Size of the internal record serialization buffer in bytes (default: 8 MB).
         *
         * @throws std::invalid_argument If input vector dimensions do not match cfg constants or buffer is too small.
         * @throws std::runtime_error If file opening, writing, or flushing fails, or if data constraints are violated.
         * @throws std::overflow_error If record counts or offsets exceed standard integer limits.
         */
        static void write_cdx_file(
            const std::filesystem::path& filepath,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            std::uint64_t nid_offset = 0,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE
        );

        /**
         * @brief Serializes graph node attributes into a Zstandard-compressed binary CDX file (.cdx.zst).
         *
         * @param filepath Target filesystem path for the compressed output file.
         * @param compo_names Vector of UTF-8 component names (must match cfg::N_COMPO).
         * @param nid2compo Mapping from node ID to component ID (size must match cfg::ARRAY_SIZE).
         * @param local_idx Mapping from node ID to component-local position (size must match cfg::ARRAY_SIZE).
         * @param seq_len Mapping from node ID to sequence length in base pairs (size must match cfg::ARRAY_SIZE).
         * @param nid_offset Offset added to array indices to compute output Node IDs (default: 0).
         * @param buffer_size Size of the internal record serialization buffer in bytes (default: 8 MB).
         * @param compression_level Zstandard compression level ranging from 1 (fastest) to 22 (highest compression,
         *          default: 3).
         *
         * @throws std::invalid_argument If compression level is outside [1, 22] or vector sizes mismatch cfg limits.
         * @throws std::runtime_error If Zstandard compression or file writing operations fail.
         * @throws std::overflow_error If uncompressed or compressed buffers exceed capacity limits.
         */
        static void write_cdx_zstd_file(
            const std::filesystem::path& filepath,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            std::uint64_t nid_offset = 0,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE,
            int compression_level = 3
        );

    private:

        /**
         * @brief Validates input array dimensions, checks sentinel consistency, and counts active records per component.
         *
         * @param compo_names UTF-8 component names.
         * @param nid2compo Component assignments per node.
         * @param local_idx Local positions per node within component.
         * @param seq_len Sequence lengths per node.
         * @return CdxWriteStats Calculated per-component record counts and total valid records.
         *
         * @throws std::invalid_argument If vector lengths do not match cfg constants or component names are invalid.
         * @throws std::runtime_error If active/inactive sentinel states are inconsistent across arrays.
         * @throws std::overflow_error If record counting overflows integer thresholds.
         */
        static CdxWriteStats validate_and_count_records(
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len
        );

#ifndef NDEBUG
        /**
         * @brief Debug-only verification step to confirm that local_idx forms a contiguous,
         * non-overlapping sequence [0, component_size).
         *
         * @note Only executed in Debug builds (compiled out when NDEBUG is defined).
         * Uses bitmap allocations to detect duplicates or missing indices.
         *
         * @param nid2compo Component assignments per node.
         * @param local_idx Local positions per node.
         * @param component_counts Total active record count per component.
         *
         * @throws std::logic_error If index duplicates, missing indices, or out-of-bounds local positions are detected.
         */
        static void debug_validate_local_idx_density(
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint64_t>& component_counts
        );
#endif

        /**
         * @brief Core low-level engine that writes the binary CDX v1 format to any output stream (file or stringstream).
         *
         * @param out Reference to an open std::ostream.
         * @param compo_names UTF-8 component names.
         * @param nid2compo Component assignments per node.
         * @param local_idx Local positions per node.
         * @param seq_len Sequence lengths per node.
         * @param nid_offset Offset added to node array index (default: 0).
         * @param buffer_size Serialization buffer capacity in bytes (default: 8 MB).
         */
        static void write_cdx_stream(
            std::ostream& out,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t>& nid2compo,
            const std::vector<std::uint32_t>& local_idx,
            const std::vector<std::uint32_t>& seq_len,
            std::uint64_t nid_offset = 0,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE
        );

        // Header tracking metrics per component (used internally during block packing).
        struct ComponentStats
        {
            std::uint64_t n_records  = 0;          ///< Total valid node records in component.
            std::uint64_t nid_min    = UINT64_MAX;  ///< Smallest node ID in component.
            std::uint64_t nid_max    = 0;           ///< Largest node ID in component.
        };
    };

    /**
     * @brief Utility class for exporting graph layout data into tab-separated values (TSV) format.
     */
    class TsvWriter
    {
    public:
        TsvWriter() = delete;
        static constexpr std::size_t DEFAULT_BUFFER_SIZE = 8 * 1024 * 1024; // 8 MB

        /**
         * @brief Writes graph node alignments and sequence coordinates to a TSV stream (e.g., std::cout or std::ofstream).
         *
         * Output Columns:
         * compo_name \t component_id \t idx \t node_id \t start_pos \t length \n
         *
         * @param out Target output stream.
         * @param compo_names Component name mapping (must match cfg::N_COMPO).
         * @param nid2compo Mapping from node ID to component ID (size must match cfg::ARRAY_SIZE).
         * @param local_idx Mapping from node ID to local linear index (size must match cfg::ARRAY_SIZE).
         * @param start_pos Mapping from node ID to calculated 0-based sequence start coordinate.
         * @param seq_len Mapping from node ID to node sequence length.
         * @param node_id_offset Offset added to array indices to compute printed Node IDs (default: 0).
         * @param buffer_size Output formatting buffer capacity in bytes (default: 8 MB).
         *
         * @throws std::invalid_argument If vector dimensions fail to match cfg constants or buffer size < 64 bytes.
         * @throws std::runtime_error If stream writing fails or invalid component identifiers are encountered.
         */
        static void write_tsv(
            std::ostream &out,
            const std::vector<std::string>& compo_names,
            const std::vector<std::uint16_t> &nid2compo,
            const std::vector<std::uint32_t> &local_idx,
            const std::vector<std::uint32_t> &start_pos,
            const std::vector<std::uint32_t> &seq_len,
            std::uint64_t node_id_offset = 0,
            std::size_t buffer_size = DEFAULT_BUFFER_SIZE
        );
    };

} // namespace cdx

#endif // CDX_BUILDER_CDX_WRITER_H