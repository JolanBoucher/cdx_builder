//
// Created by Jolan on 2026-07-23.
//

#ifndef CDX_BUILDER_CONSTANTS_H
#define CDX_BUILDER_CONSTANTS_H

#include <cstdint>

namespace cfg {

    // Sentinels using to tell that a nid don't represent any node
    inline constexpr uint32_t NODE_UNSEEN_32 = 0xFFFFFFFF;
    inline constexpr uint16_t NODE_UNSEEN_16 = 0xFFFF;
    inline constexpr uint8_t NODE_UNSEEN_8 = 0xFF;

    // Values initialized after loading the graph.
    inline size_t NB_NODES = 0;
    inline size_t ARRAY_SIZE = 0;
    inline size_t N_COMPO = 0;
    inline size_t N_HAPLO = 0;


}

#endif