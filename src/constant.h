/**
 * @file constant.h
 * @brief Global configuration constants and runtime parameters for the CDX builder.
 *
 * This module defines sentinel values for uninitialized or missing node identifiers
 * across different integer widths, alongside mutable global state variables initialized
 * during graph loading (such as total node counts, array sizing bounds, and haplotype statistics).
 */

#ifndef CDX_BUILDER_CONSTANTS_H
#define CDX_BUILDER_CONSTANTS_H

#include <cstdint>
#include <cstddef>

namespace cfg {

    // --- Sentinel Constants ---
    inline constexpr uint32_t NODE_UNSEEN_32 = 0xFFFFFFFF;  // Sentinel representing an unused or invalid 32-bit node ID
    inline constexpr uint16_t NODE_UNSEEN_16 = 0xFFFF;      // Sentinel representing an unused or invalid 16-bit node ID
    inline constexpr uint8_t  NODE_UNSEEN_8  = 0xFF;        // Sentinel representing an unused or invalid 8-bit node ID

    // --- Runtime Graph Parameters (Initialized Post-Loading) ---
    inline size_t NB_NODES    = 0;      // Total number of active nodes loaded in the graph
    inline size_t ARRAY_SIZE  = 0;      // Allocated capacity for node-indexed global vectors
    inline size_t N_COMPO     = 0;      // Total number of connected components
    inline size_t N_HAPLO     = 0;      // Total number of haplotypes in the graph

} // namespace cfg

#endif // CDX_BUILDER_CONSTANTS_H