# ------------------------------------------------------------------------------
# Compiler requirements
# ------------------------------------------------------------------------------
#
# src/ relies on C++17 language/library features that need a reasonably recent
# GCC. Apple Clang and other compilers are left unchecked here (Apple Clang's
# version numbering doesn't map to upstream Clang releases, and it already
# supports everything we need on any macOS version CMake itself supports).

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
    message(FATAL_ERROR
            "GCC/G++ 11 or newer is required. "
            "Detected version: ${CMAKE_CXX_COMPILER_VERSION}")
endif()
