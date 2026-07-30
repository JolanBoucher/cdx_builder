//
// Created by Jolan on 2026-07-27.
//

#include "cdx_format.h"

namespace cdx
{
    // Converts a 32-bit unsigned integer to little-endian byte order.
    uint32_t CdxFormat::to_little_endian32(const uint32_t value) noexcept {
        if constexpr (IS_LITTLE_ENDIAN) return value;

    #if defined(__GNUC__) || defined(__clang__)
        return __builtin_bswap32(value);

    #else
        return
            ((value & 0x000000FFU) << 24U) |
            ((value & 0x0000FF00U) << 8U)  |
            ((value & 0x00FF0000U) >> 8U)  |
            ((value & 0xFF000000U) >> 24U);
    #endif
    }

    // Converts a 64-bit unsigned integer to little-endian byte order.
    uint64_t CdxFormat::to_little_endian64(const uint64_t value) noexcept {
        if constexpr (IS_LITTLE_ENDIAN) return value;
    #if defined(__GNUC__) || defined(__clang__)
        return __builtin_bswap64(value);
    #else
        return
            ((value & 0x00000000000000FFULL) << 56U) |
            ((value & 0x000000000000FF00ULL) << 40U) |
            ((value & 0x0000000000FF0000ULL) << 24U) |
            ((value & 0x00000000FF000000ULL) << 8U)  |
            ((value & 0x000000FF00000000ULL) >> 8U)  |
            ((value & 0x0000FF0000000000ULL) >> 24U) |
            ((value & 0x00FF000000000000ULL) >> 40U) |
            ((value & 0xFF00000000000000ULL) >> 56U);
    #endif
    }

    // Converts a 32-bit little-endian integer to native byte order.
    uint32_t CdxFormat::from_little_endian32(const uint32_t value) noexcept{
        return to_little_endian32(value);
    }

    // Converts a 64-bit little-endian integer to native byte order.
    uint64_t CdxFormat::from_little_endian64(const uint64_t value) noexcept{
        return to_little_endian64(value);
    }

    // Serializes the global file header into the destination buffer with little-endian encoding.
    void CdxFormat::pack_file_header(
        char *destination,
        const uint32_t component_count) noexcept {
        CdxFileHeader header{};
        header.magic        = MAGIC;
        header.n_components = to_little_endian32(component_count);
        header.nid_width    = NID_WIDTH;
        header.seqlen_width = SEQLEN_WIDTH;

        std::memcpy(destination, &header, sizeof(header));
    }

    // Serializes a component block header into the destination buffer with little-endian encoding.
    void CdxFormat::pack_component_header(
        char* destination,
        const uint64_t record_count,
        const uint64_t minimum_node_id,
        const uint64_t maximum_node_id,
        const uint32_t name_size) noexcept {
        CdxComponentHeader header{};
        header.n_records   = to_little_endian64(record_count);
        header.node_id_min = to_little_endian64(minimum_node_id);
        header.node_id_max = to_little_endian64(maximum_node_id);
        header.name_size   = to_little_endian32(name_size);

        std::memcpy(destination, &header, sizeof(header));
    }

    // Serializes an individual node record into the destination buffer with little-endian encoding.
    void CdxFormat::pack_node_record(
        char *destination,
        const uint64_t node_id,
        const uint32_t local_index,
        const uint32_t sequence_length) noexcept {
        CdxNodeRecord record{};
        record.node_id  = to_little_endian64(node_id);
        record.idx      = to_little_endian32(local_index);
        record.seq_len  = to_little_endian32(sequence_length);

        std::memcpy(destination, &record, sizeof(record));
    }

    // Deserializes and converts the global file header from the source buffer into native byte order.
    CdxFileHeader CdxFormat::unpack_file_header(const char* source) noexcept {
        CdxFileHeader header{};
        std::memcpy(&header, source, sizeof(header));
        header.n_components  = from_little_endian32(header.n_components);
        return header;
    }

    // Deserializes and converts a component header from the source buffer into native byte order.
    CdxComponentHeader CdxFormat::unpack_component_header(const char* source) noexcept {
        CdxComponentHeader header{};
        std::memcpy(&header, source, sizeof(header));
        header.n_records    = from_little_endian64(header.n_records);
        header.node_id_min  = from_little_endian64(header.node_id_min);
        header.node_id_max  = from_little_endian64(header.node_id_max);
        header.name_size    = from_little_endian32(header.name_size);
        return header;
    }

    // Deserializes and converts a node record from the source buffer into native byte order.
    CdxNodeRecord CdxFormat::unpack_node_record(const char* source) noexcept {
        CdxNodeRecord record{};
        std::memcpy(&record, source, sizeof(record));
        record.node_id      = from_little_endian64(record.node_id);
        record.idx          = from_little_endian32(record.idx);
        record.seq_len      = from_little_endian32(record.seq_len);
        return record;
    }

    // Validates whether the file header contains the expected CDX magic signature bytes.
    bool CdxFormat::has_valid_magic(const CdxFileHeader& header) noexcept {
        return header.magic == MAGIC;
    }

    // Validates whether the file header specifies supported node ID and sequence length byte widths.
    bool CdxFormat::has_valid_widths(const CdxFileHeader& header) noexcept {
        return header.nid_width == NID_WIDTH &&
               header.seqlen_width == SEQLEN_WIDTH;
    }
}