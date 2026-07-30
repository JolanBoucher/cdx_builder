/**
 * @file cdx_format.h
 * @brief CDX Archive Format Specification and Binary Serialization Module.
 *
 * This header defines the complete binary layout, validation assertions, and serialization
 * specification for CDX archive files. It provides:
 *
 * 1. **Packed Binary Record Structures**:
 *    - Strict 1-byte aligned structs (`CdxFileHeader`, `CdxComponentHeader`, `CdxNodeRecord`)
 *      guaranteed by compile-time assertions to match exact disk footprints.
 *
 * 2. **`CdxFormat` Static Utility Class**:
 *    - Centralizes archive type aliases, magic byte signatures (`MAGIC`), fixed byte-width specs,
 *      and component structure sizing parameters.
 *    - Enforces cross-platform compatibility via compile-time endian detection (`IS_LITTLE_ENDIAN`)
 *      and optimized byte-swapping utility routines.
 *    - Exposes zero-copy buffer serialization (`pack_*`) and deserialization (`unpack_*`) APIs
 *      alongside format integrity validation helpers.
 */

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
     * @brief Global 10-byte archive header.
     *
     * Located at the very beginning of a CDX file to validate file format integrity
     * and specify global metadata dimensions.
     *
     * Layout: magic (4B) | n_components (4B) | nid_width (1B) | seqlen_width (1B)
     */
    struct CdxFileHeader
    {
        std::array<char, 4> magic;      // Magic number/bytes identifying the valid CDX file format signature
        uint32_t n_components;          // Total number of connected components stored in the archive
        uint8_t nid_width;              // Byte width descriptor for node identifiers (e.g., 1, 2, 4, or 8 bytes)
        uint8_t seqlen_width;           // Byte width descriptor for sequence length fields
    };


    /**
     * @brief Localized 28-byte block header for an individual graph component.
     *
     * Precedes each component's data section within the archive to describe
     * the range, size, and metadata name length of the contained nodes.
     *
     * Layout: n_records (8B) | node_id_min (8B) | node_id_max (8B) | name_size (4B)
     */
    struct CdxComponentHeader
    {
        uint64_t n_records;             // Total number of node records belonging to this component
        uint64_t node_id_min;           // Minimum node identifier value present in this component block
        uint64_t node_id_max;           // Maximum node identifier value present in this component block
        uint32_t name_size;             // Byte size of the component's name or string identifier string
    };

    /**
     * @brief Continuous 16-byte topological node record.
     *
     * Represents the serialized spatial and indexing attributes of a single graph node.
     *
     * Layout: node_id (8B) | idx (4B) | seq_len (4B)
     */
    struct CdxNodeRecord
    {
        uint64_t node_id;           // Unique global identifier of the graph node
        uint32_t idx;               // Assigned local sequential rank or layout index within the component
        uint32_t seq_len;           // Biological sequence length associated with the node
    };
    #pragma  pack(pop)

    // Static assertion checks to guarantee exact binary layout size
    static_assert(sizeof(CdxFileHeader) == 10, "CdxFileHeader must occupy exactly 10 bytes");
    static_assert(sizeof(CdxComponentHeader) == 28, "CdxComponentHeader must occupy exactly 28 bytes");
    static_assert(sizeof(CdxNodeRecord) == 16, "CdxNodeRecord must occupy exactly 16 bytes");

    static_assert(std::is_trivially_copyable_v<CdxFileHeader>, "CdxFileHeader must be trivially copyable");
    static_assert(std::is_trivially_copyable_v<CdxComponentHeader>, "CdxComponentHeader must be trivially copyable");
    static_assert(std::is_trivially_copyable_v<CdxNodeRecord>, "CdxNodeRecord must be trivially copyable");

    /**
     * @brief Encapsulates binary format specifications, type aliases, serialization routines,
     *        and endianness utilities for CDX graph archive generation and parsing.
     *
     * `CdxFormat` acts as a static utility class that defines layout constants, block sizes,
     * and compile-time endian detection, providing high-performance methods to pack and unpack
     * file headers, component headers, and node records directly from raw byte buffers.
     */
    class CdxFormat final
    {
    public:
        // --- Type Aliases ---
        using ComponentId   = std::uint64_t;        // Type alias for connected component identifiers
        using Nid           = std::uint64_t;        // Type alias for global node identifiers
        using Idx           = std::uint32_t;        // Type alias for localized local layout indices
        using SeqLen        = std::uint32_t;        // Type alias for biological sequence lengths
        using RecordCount   = std::uint64_t;        // Type alias for record count counters

        // --- Format Constants ---
        static constexpr std::array<char, 4> MAGIC = {'C', 'D', 'X', '\x01'};               // Magic signature bytes identifying valid CDX files
        static constexpr uint8_t NID_WIDTH  = 8;                                            // Standard byte width for serialized node identifiers (8 bytes)
        static constexpr uint8_t SEQLEN_WIDTH  = 4;                                         // Standard byte width for serialized sequence lengths (4 bytes)

        // --- Serialized Block Size Descriptors ---
        static constexpr std::size_t FILE_HEADER_SIZE       = sizeof(CdxFileHeader);       // Total byte size of the global file header (10 bytes)
        static constexpr std::size_t COMPONENT_HEADER_SIZE  = sizeof(CdxComponentHeader);  // Total byte size of each component block header (28 bytes)
        static constexpr std::size_t RECORD_SIZE            = sizeof(CdxNodeRecord);       // Total byte size of an individual node record (16 bytes)

        /**
         * @brief Strict compile-time endianness detection block.
         *
         * Evaluates standard compiler preprocessor macros, target operating systems,
         * and explicit hardware architecture flags to determine whether the host system
         * uses little-endian byte ordering (`IS_LITTLE_ENDIAN = true`) or big-endian
         * byte ordering (`IS_LITTLE_ENDIAN = false`). Halts compilation if the byte order
         * cannot be reliably resolved.
         */
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

        // --- Constructor Deletion ---
        CdxFormat() = delete;


        // --- Endianness Conversion Utilities ---
        [[nodiscard]] static uint32_t to_little_endian32(uint32_t value) noexcept;
        [[nodiscard]] static uint64_t to_little_endian64(uint64_t value) noexcept;
        [[nodiscard]] static uint32_t from_little_endian32(uint32_t value) noexcept;
        [[nodiscard]] static uint64_t from_little_endian64(uint64_t value) noexcept;

        // --------------------------------------- //
        // --- Serialization (Packing) Methods --- //
        // --------------------------------------- //
        /**
         * @brief Serializes the global file header into the destination buffer in little-endian format.
         * @pre destination points to at least FILE_HEADER_SIZE writable bytes.
         */
        static void pack_file_header(
            char* destination,
            uint32_t component_count) noexcept;

        /**
         * @brief Serializes a component block header into the destination buffer in little-endian format.
         * @pre destination points to at least COMPONENT_HEADER_SIZE writable bytes.
         */
        static void pack_component_header(
            char* destination,
            uint64_t record_count,
            uint64_t minimum_node_id,
            uint64_t maximum_node_id,
            uint32_t name_size) noexcept;

        /**
         * @brief Serializes an individual node record into the destination buffer in little-endian format.
         * @pre destination points to at least RECORD_SIZE writable bytes.
         */
        static void pack_node_record(
            char* destination,
            uint64_t node_id,
            uint32_t local_index,
            uint32_t sequence_length) noexcept;

        // ------------------------------------------- //
        // --- Deserialization (Unpacking) Methods --- //
        // ------------------------------------------- //
        /**
         * @brief Unpacks and converts a global file header from the source buffer, adjusting endianness.
         * @pre source points to at least FILE_HEADER_SIZE readable bytes.
         */
        [[nodiscard]] static CdxFileHeader unpack_file_header(const char* source) noexcept;

        /**
         * @brief Unpacks and converts a component header from the source buffer, adjusting endianness.
         * @pre source points to at least COMPONENT_HEADER_SIZE readable bytes.
         */
        [[nodiscard]] static CdxComponentHeader unpack_component_header(const char* source) noexcept;

        /**
         * @brief Unpacks and converts a node record from the source buffer, adjusting endianness.
         * @pre source points to at least RECORD_SIZE readable bytes.
         */
        [[nodiscard]] static CdxNodeRecord unpack_node_record(const char* source) noexcept;

        // --- Validation Helpers ---
        [[nodiscard]] static bool has_valid_magic(const CdxFileHeader& header) noexcept;
        [[nodiscard]] static bool has_valid_widths(const CdxFileHeader& header) noexcept;
    };
} // namespace cdx
#endif //CDX_BUILDER_CDX_FORMAT_H