//
// Created by Jolan on 2026-07-27.
//
#ifndef CDX_BUILDER_CDX_FORMAT_H
#define CDX_BUILDER_CDX_FORMAT_H

#include <array>
#include <cstddef>
#include <cstring>
#include <cstdint>
#include <type_traits>

namespace cdx {
#pragma pack(push, 1)
    /**
    * Global 10-byte archive header.
    * Layout: magic (4B) | n_components (4B) | nid_width (1B) | seqlen_width (1B)
    */
    struct CdxFileHeader
    {
        std::array<char, 4> magic;
        uint32_t n_components;
        uint8_t nid_width;
        uint8_t seqlen_width;
    };

    /**
     * Localized 28-byte block header for a component.
     * Layout: n_records (8B) | node_id_min (8B) | node_id_max (8B) | name_size (4B)
     */
    struct CdxComponentHeader
    {
        uint64_t n_records;
        uint64_t node_id_min;
        uint64_t node_id_max;
        uint32_t name_size;
    };

    /**
     * Continuous 16-byte topological node record.
     * Layout: node_id (8B) | idx (4B) | seq_len (4B)
     */
    struct CdxNodeRecord
    {
        uint64_t node_id;
        uint32_t idx;
        uint32_t seq_len;
    };

    #pragma  pack(pop)

    // Static assertion checks to guarantee exact binary layout size
    static_assert(sizeof(CdxFileHeader) == 10, "CdxFileHeader must occupy exactly 10 bytes");
    static_assert(sizeof(CdxComponentHeader) == 28, "CdxComponentHeader must occupy exactly 28 bytes");
    static_assert(sizeof(CdxNodeRecord) == 16, "CdxNodeRecord must occupy exactly 16 bytes");

    static_assert(std::is_trivially_copyable_v<CdxFileHeader>, "CdxFileHeader must be trivially copyable");
    static_assert(std::is_trivially_copyable_v<CdxComponentHeader>, "CdxComponentHeader must be trivially copyable");
    static_assert(std::is_trivially_copyable_v<CdxNodeRecord>, "CdxNodeRecord must be trivially copyable");

    class CdxFormat final
    {
    public:
        // constants
        using ComponentId   = std::uint64_t;
        using Nid           = std::uint64_t;
        using Idx           = std::uint32_t;
        using SeqLen        = std::uint32_t;
        using RecordCount   = std::uint64_t;

        static constexpr std::array<char, 4> MAGIC = {'C', 'D', 'X', '\x01'};
        static constexpr uint8_t NID_WIDTH  = 8;
        static constexpr uint8_t SEQLEN_WIDTH  = 4;

        static constexpr std::size_t FILE_HEADER_SIZE  = sizeof(CdxFileHeader);
        static constexpr std::size_t COMPONENT_HEADER_SIZE  = sizeof(CdxComponentHeader);
        static constexpr std::size_t RECORD_SIZE  = sizeof(CdxNodeRecord);


        // Endian detection
        // Strict compile-time detection of the native byte order.
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && defined(__ORDER_BIG_ENDIAN__)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    inline static constexpr bool IS_LITTLE_ENDIAN = true;
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    inline static constexpr bool IS_LITTLE_ENDIAN = false;
#else
#error "Unsupported mixed-endian architecture"
#endif
#elif defined(_WIN32)
    inline static constexpr bool IS_LITTLE_ENDIAN = true;
#elif defined(__x86_64__) || defined(__i386__) || defined(__AARCH64EL__) || defined(__ARMEL__)
    inline static constexpr bool IS_LITTLE_ENDIAN = true;     // Common explicitly little-endian architectures.
#elif defined(__AARCH64EB__) || defined(__ARMEB__) || defined(__s390x__)
    inline static constexpr bool IS_LITTLE_ENDIAN = false;    // Common explicitly big-endian architectures.
#else
#error "Unable to determine native byte order for the CDX format"
#endif

        // no constructor
        CdxFormat() = delete;


        [[nodiscard]] static uint32_t to_little_endian32(uint32_t value) noexcept;
        [[nodiscard]] static uint64_t to_little_endian64(uint64_t value) noexcept;
        [[nodiscard]] static uint32_t from_little_endian32(uint32_t value) noexcept;
        [[nodiscard]] static uint64_t from_little_endian64(uint64_t value) noexcept;

        // Methods pack/unpack
        /**
         * Precondition: destination points to at least FILE_HEADER_SIZE writable bytes.
         */
        static void pack_file_header(
            char* destination,
            uint32_t component_count) noexcept;

        /**
         * Precondition: destination points to at least COMPONENT_HEADER_SIZE writable bytes.
         */
        static void pack_component_header(
            char* destination,
            uint64_t record_count,
            uint64_t minimum_node_id,
            uint64_t maximum_node_id,
            uint32_t name_size) noexcept;

        /**
         * Precondition: destination points to at least RECORD_SIZE writable bytes.
         */
        static void pack_node_record(
            char* destination,
            uint64_t node_id,
            uint32_t local_index,
            uint32_t sequence_length) noexcept;

        /**
         * Precondition: source points to at least FILE_HEADER_SIZE readable bytes.
         */
        [[nodiscard]] static CdxFileHeader unpack_file_header(const char* source) noexcept;

        /**
         * Precondition: source points to at least COMPONENT_HEADER_SIZE readable bytes.
         */
        [[nodiscard]] static CdxComponentHeader unpack_component_header(const char* source) noexcept;

        /**
         * Precondition: source points to at least RECORD_SIZE readable bytes.
         */
        [[nodiscard]] static CdxNodeRecord unpack_node_record(const char* source) noexcept;

        [[nodiscard]] static bool has_valid_magic(const CdxFileHeader& header) noexcept;
        [[nodiscard]] static bool has_valid_widths(const CdxFileHeader& header) noexcept;
    };
} // namespace cdx
#endif //CDX_BUILDER_CDX_FORMAT_H