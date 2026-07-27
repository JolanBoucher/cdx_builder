cmake_minimum_required(VERSION 3.16)
project(cdx_builder CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# 0. System Dependencies (OpenMP, Threads, OpenSSL, zstd)
# ------------------------------------------------------------------------------
if(APPLE)
    # macOS-specific Homebrew paths
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")

    list(APPEND CMAKE_PREFIX_PATH
            "/opt/homebrew"
            "/opt/homebrew/opt/openssl@3"
            "/opt/homebrew/opt/zstd"
    )
endif()

# Standard dependency detection (Works natively on Linux!)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)

find_library(ZSTD_LIBRARY NAMES zstd REQUIRED)
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h REQUIRED)

# ------------------------------------------------------------------------------
# 1. External Dependency: libbdsg
# ------------------------------------------------------------------------------
set(BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "" FORCE)

add_subdirectory(deps/libbdsg)

# ------------------------------------------------------------------------------
# 2. External Dependencies: gbwt & gbwtgraph
# ------------------------------------------------------------------------------
include(ExternalProject)

set(GBWT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/gbwt")
set(GBWTGRAPH_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/gbwtgraph")
set(LIBBDSG_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps/libbdsg")

set(HANDLEGRAPH_INCLUDE_DIR "${LIBBDSG_DIR}/bdsg/deps/libhandlegraph/src/include")
set(HANDLEGRAPH_BUILD_DIR "${LIBBDSG_DIR}/bdsg/deps/libhandlegraph/build")

set(GBWT_LIBRARY "${GBWT_DIR}/lib/libgbwt.a")
set(GBWTGRAPH_LIBRARY "${GBWTGRAPH_DIR}/lib/libgbwtgraph.a")

set(GBWT_INCLUDE_DIR "${GBWT_DIR}/include")
set(GBWTGRAPH_INCLUDE_DIR "${GBWTGRAPH_DIR}/include")

ExternalProject_Add(
        gbwt_ext
        SOURCE_DIR "${GBWT_DIR}"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -j4
        INSTALL_COMMAND ""
        BUILD_BYPRODUCTS "${GBWT_LIBRARY}"
)

ExternalProject_Add(
        gbwtgraph_ext
        DEPENDS gbwt_ext
        SOURCE_DIR "${GBWTGRAPH_DIR}"
        BUILD_IN_SOURCE 1
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -j4 GBWT_DIR=${GBWT_DIR} CPPFLAGS=-I${HANDLEGRAPH_INCLUDE_DIR} LDFLAGS=-L${HANDLEGRAPH_BUILD_DIR}
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

# ------------------------------------------------------------------------------
# 3. Main Executable Target
# ------------------------------------------------------------------------------
add_executable(cdx_builder
        src/main.cpp
        src/gbz_IO.cpp
        src/gbz_IO.h
        src/constant.h
        src/graph_linearization.cpp
        src/graph_linearization.h
        src/CSRMatrix.cpp
        src/graph_indexing.cpp
        src/graph_indexing.h
        src/cdx_format.cpp
        src/cdx_format.h
)

target_include_directories(cdx_builder PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${GBWT_INCLUDE_DIR}
        ${GBWTGRAPH_INCLUDE_DIR}
        ${ZSTD_INCLUDE_DIR}
)

target_link_libraries(cdx_builder PRIVATE
        libbdsg
        gbwtgraph_lib
        gbwt_lib
        OpenSSL::Crypto
        ${ZSTD_LIBRARY}
        OpenMP::OpenMP_CXX
        Threads::Threads
        z
)

# ------------------------------------------------------------------------------
# 4. Google Test
# ------------------------------------------------------------------------------
enable_testing()

add_executable(run_tests
        tests/test_gbz_io.cpp
        src/gbz_IO.cpp
)

target_include_directories(run_tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${GBWT_INCLUDE_DIR}
        ${GBWTGRAPH_INCLUDE_DIR}
        ${ZSTD_INCLUDE_DIR}
)

target_link_libraries(run_tests PRIVATE
        gtest_main
        gtest
        libbdsg
        gbwtgraph_lib
        gbwt_lib
        OpenSSL::Crypto
        ${ZSTD_LIBRARY}
        OpenMP::OpenMP_CXX
        Threads::Threads
        z
)

include(GoogleTest)
gtest_discover_tests(run_tests)