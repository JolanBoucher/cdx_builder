# ------------------------------------------------------------------------------
# External dependencies
# ------------------------------------------------------------------------------
#
# Resolves every third-party dependency needed by the cdx_builder target:
#   1. System libraries expected to already be installed (OpenMP, OpenSSL,
#      zstd, Boost, Jansson, CLI11).
#   2. libbdsg, vendored as a submodule and built directly via add_subdirectory.
#   3. sdsl-lite, GBWT and GBWTGraph, built out-of-tree via ExternalProject
#      (they don't build cleanly as plain add_subdirectory targets alongside
#      libbdsg's own copy of sdsl/handlegraph).
#   4. cdx_lib, the CDX format/IO code shared with cdx_coverage.
#
# All variables set here (GBWT_INCLUDE_DIR, GBWTGRAPH_INCLUDE_DIR,
# ZSTD_LIBRARY, ZSTD_INCLUDE_DIR, JANSSON_*, LIBBDSG_DIR, imported targets
# gbwt_lib/gbwtgraph_lib/libbdsg/cdx_lib, ...) are consumed by the root
# CMakeLists.txt and by tests/CMakeLists.txt through normal directory-scope
# inheritance.

# --- 1. System libraries ------------------------------------------------------

find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Boost REQUIRED)

find_library(ZSTD_LIBRARY NAMES zstd REQUIRED)
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h REQUIRED)

pkg_check_modules(JANSSON REQUIRED jansson)

message(STATUS "Found OpenMP")
message(STATUS "Found OpenSSL")
message(STATUS "Found zstd")
message(STATUS "Found Boost")
message(STATUS "Found Jansson")

add_subdirectory(deps/CLI11)

# --- 2. libbdsg ---------------------------------------------------------------

# Disable tests and optional targets from bundled dependencies.
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(SDSL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "" FORCE)

add_subdirectory(deps/libbdsg)

# --- 3. sdsl-lite, GBWT, GBWTGraph --------------------------------------------

include(ExternalProject)

set(GBWT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/gbwt")
set(GBWTGRAPH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/gbwtgraph")
set(LIBBDSG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libbdsg")
set(SDSL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/sdsl-lite")

ExternalProject_Add(
        sdsl_ext
        SOURCE_DIR "${SDSL_DIR}"
        BINARY_DIR "${CMAKE_BINARY_DIR}/sdsl_ext-build"
        CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
        -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/sdsl-install
        -DCMAKE_POLICY_VERSION_MINIMUM=3.5
        INSTALL_COMMAND ${CMAKE_COMMAND} --build . --target install -- -j1
)

set(HANDLEGRAPH_INCLUDE_DIR
        "${LIBBDSG_DIR}/bdsg/deps/libhandlegraph/src/include")
set(HANDLEGRAPH_BUILD_DIR
        "${LIBBDSG_DIR}/bdsg/deps/libhandlegraph/build")

set(GBWT_LIBRARY "${GBWT_DIR}/lib/libgbwt.a")
set(GBWTGRAPH_LIBRARY "${GBWTGRAPH_DIR}/lib/libgbwtgraph.a")
set(GBWT_INCLUDE_DIR "${GBWT_DIR}/include")
set(GBWTGRAPH_INCLUDE_DIR "${GBWTGRAPH_DIR}/include")

# ---- GBWT ------------------------------------------------------

ExternalProject_Add(
        gbwt_ext
        DEPENDS sdsl_ext
        SOURCE_DIR "${GBWT_DIR}"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
        ${CMAKE_COMMAND} -E env
        CPATH=${CMAKE_CURRENT_SOURCE_DIR}/deps/sdsl-lite/include
        LIBRARY_PATH=${LIBBDSG_DIR}/lib
        LD_LIBRARY_PATH=${LIBBDSG_DIR}/lib
        make -j4
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${GBWT_LIBRARY}"
)

# ---- GBWTGRAPH -------------------------------------------------

ExternalProject_Add(
        gbwtgraph_ext
        DEPENDS gbwt_ext
        SOURCE_DIR "${GBWTGRAPH_DIR}"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND
        ${CMAKE_COMMAND} -E env
        CPATH=${GBWT_INCLUDE_DIR}:${HANDLEGRAPH_INCLUDE_DIR}:/usr/local/include
        LIBRARY_PATH=${LIBBDSG_DIR}/lib:${GBWT_DIR}/lib:/usr/local/lib
        make -j4
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${GBWTGRAPH_LIBRARY}"
)

add_dependencies(gbwtgraph_ext libbdsg)

add_library(gbwt_lib STATIC IMPORTED)
set_target_properties(gbwt_lib PROPERTIES IMPORTED_LOCATION "${GBWT_LIBRARY}")
add_dependencies(gbwt_lib gbwt_ext)

add_library(gbwtgraph_lib STATIC IMPORTED)
set_target_properties(gbwtgraph_lib PROPERTIES IMPORTED_LOCATION "${GBWTGRAPH_LIBRARY}")
add_dependencies(gbwtgraph_lib gbwtgraph_ext)

# --- 4. cdx_lib -----------------------------------------------------------
#
# cdx_lib is shared between cdx_builder and cdx_coverage and is vendored here
# as a git submodule under deps/cdx_lib (same convention as the other deps).
#
# If you cloned this repository without submodules, run:
#   git submodule update --init --recursive
#
# CDX_LIB_DIR can be overridden (e.g. -DCDX_LIB_DIR=/path/to/cdx_lib) for
# local development against a working copy that lives elsewhere.
set(CDX_LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/cdx_lib" CACHE PATH
        "Path to the cdx_lib source tree")

if(NOT EXISTS "${CDX_LIB_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
            "cdx_lib not found at ${CDX_LIB_DIR}.\n"
            "Fetch it with:\n"
            "    git submodule update --init --recursive\n"
            "or point CDX_LIB_DIR at a local checkout, e.g.:\n"
            "    cmake -S . -B build -DCDX_LIB_DIR=/path/to/cdx_lib")
endif()

add_subdirectory("${CDX_LIB_DIR}" cdx_lib_build)
