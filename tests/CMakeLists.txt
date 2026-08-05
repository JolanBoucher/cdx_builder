# ------------------------------------------------------------------------------
# Unit test suites for cdx_builder
# ------------------------------------------------------------------------------
#
# Every module under src/ that can be exercised in isolation gets its own test
# executable here, linked only against the dependencies it actually needs
# (e.g. test_csr_matrix does not need to link gbwt/gbwtgraph/libbdsg at all).
# This keeps incremental rebuilds fast and keeps a ctest failure scoped to a
# single module instead of one monolithic test binary.
#
# All variables referenced below (GBWT_INCLUDE_DIR, ZSTD_LIBRARY, imported
# targets such as gbwt_lib/gbwtgraph_lib/libbdsg/cdx_lib, etc.) are defined in
# the parent CMakeLists.txt and are visible here through normal CMake
# directory-scope inheritance.

include(GoogleTest)

# add_cdx_test(<name>
#     SOURCES  <gtest source file(s) + the src/*.cpp file(s) under test>
#     INCLUDES <extra include directories, optional>
#     LIBS     <extra libraries to link, optional>
# )
#
# Always links against gtest/gtest_main and against ${CMAKE_SOURCE_DIR}/src so
# test files can '#include "some_module.h"' directly.
function(add_cdx_test NAME)
    cmake_parse_arguments(ARG "" "" "SOURCES;INCLUDES;LIBS" ${ARGN})

    # The gtest driver file is always the first SOURCES entry. Skip the suite
    # gracefully if it hasn't been written yet instead of failing CMake
    # generation for the whole project.
    list(GET ARG_SOURCES 0 ARG_DRIVER)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${ARG_DRIVER}")
        message(STATUS "Skipping test '${NAME}': ${ARG_DRIVER} not found yet")
        return()
    endif()

    add_executable(${NAME} ${ARG_SOURCES})

    target_include_directories(${NAME} PRIVATE
            ${CMAKE_SOURCE_DIR}/src
            ${ARG_INCLUDES}
    )

    target_link_libraries(${NAME} PRIVATE
            gtest_main
            gtest
            ${ARG_LIBS}
    )

    # Needed by any test that (transitively) links libbdsg's shared objects.
    set_target_properties(${NAME} PROPERTIES
            BUILD_RPATH "${LIBBDSG_DIR}/lib"
    )

    gtest_discover_tests(${NAME})
endfunction()

# --- CSRMatrix: pure container/invariants, no external dependencies ----------
add_cdx_test(test_csr_matrix
        SOURCES
        test_csr_matrix.cpp
        ${CMAKE_SOURCE_DIR}/src/CSRMatrix.cpp
)

# --- graph_linearization: CSR construction + topological relaxation ----------
add_cdx_test(test_graph_linearization
        SOURCES
        test_graph_linearization.cpp
        ${CMAKE_SOURCE_DIR}/src/graph_linearization.cpp
        ${CMAKE_SOURCE_DIR}/src/CSRMatrix.cpp
        INCLUDES
        ${Boost_INCLUDE_DIRS}
        LIBS
        OpenMP::OpenMP_CXX
        Threads::Threads
)

# --- graph_indexing: midpoint calculation + local rank assignment ------------
add_cdx_test(test_graph_indexing
        SOURCES
        test_graph_indexing.cpp
        ${CMAKE_SOURCE_DIR}/src/graph_indexing.cpp
)

# --- cdx_writer: binary CDX + TSV serialization -------------------------------
add_cdx_test(test_cdx_writer
        SOURCES
        test_cdx_writer.cpp
        ${CMAKE_SOURCE_DIR}/src/cdx_writer.cpp
        INCLUDES
        ${ZSTD_INCLUDE_DIR}
        LIBS
        cdx_lib
        ${ZSTD_LIBRARY}
)

# --- cli: argument parsing and output path resolution -------------------------
add_cdx_test(test_cli
        SOURCES
        test_cli.cpp
        ${CMAKE_SOURCE_DIR}/src/cli.cpp
        LIBS
        CLI11::CLI11
)

# --- gbz_IO: GBZ/GBWTGraph loading, component & path utilities ---------------
add_cdx_test(test_gbz_io
        SOURCES
        test_gbz_io.cpp
        ${CMAKE_SOURCE_DIR}/src/gbz_IO.cpp
        INCLUDES
        ${GBWT_INCLUDE_DIR}
        ${GBWTGRAPH_INCLUDE_DIR}
        ${ZSTD_INCLUDE_DIR}
        LIBS
        libbdsg
        gbwtgraph_lib
        gbwt_lib
        OpenSSL::Crypto
        ${ZSTD_LIBRARY}
        OpenMP::OpenMP_CXX
        Threads::Threads
        z
)
