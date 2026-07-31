cmake_minimum_required(VERSION 3.16)
project(cdx_builder CXX)

# ------------------------------------------------------------------------------
# 0. Build Configuration
# ------------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ------------------------------------------------------------------------------
# 1. Compiler Requirements
# ------------------------------------------------------------------------------

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
        message(FATAL_ERROR
                "GCC/G++ 11 or newer is required. "
                "Detected version: ${CMAKE_CXX_COMPILER_VERSION}")
    endif()
endif()

# ------------------------------------------------------------------------------
# 2. System Dependencies
# ------------------------------------------------------------------------------

if(APPLE)
    # macOS-specific Homebrew paths
    set(OpenMP_C_FLAGS
            "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")
    set(OpenMP_CXX_FLAGS
            "-Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include")

    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")

    set(OpenMP_omp_LIBRARY
            "/opt/homebrew/opt/libomp/lib/libomp.dylib")

    list(APPEND CMAKE_PREFIX_PATH
            "/opt/homebrew"
            "/opt/homebrew/opt/openssl@3"
            "/opt/homebrew/opt/zstd"
            "/opt/homebrew/opt/cli11"
            "/opt/homebrew/opt/jansson"
    )

    # Assure que le linker trouve les libs installées par brew dans /opt/homebrew/lib
    link_directories(/opt/homebrew/lib)
endif()

# Standard dependencies
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)
find_package(OpenSSL REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Boost REQUIRED)

# zstd
find_library(ZSTD_LIBRARY NAMES zstd REQUIRED)
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h REQUIRED)

# jansson
pkg_check_modules(JANSSON REQUIRED jansson)

message(STATUS "Found OpenMP")
message(STATUS "Found OpenSSL")
message(STATUS "Found zstd")
message(STATUS "Found Boost")
message(STATUS "Found Jansson")

# CLI11
add_subdirectory(deps/CLI11)

# ------------------------------------------------------------------------------
# 3. External Dependency: libbdsg
# ------------------------------------------------------------------------------

# Disable tests and optional targets from bundled dependencies.
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
set(SDSL_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON OFF CACHE BOOL "" FORCE)
set(BUILD_PYTHON_BINDINGS OFF CACHE BOOL "" FORCE)
set(RUN_DOXYGEN OFF CACHE BOOL "" FORCE)

add_subdirectory(deps/libbdsg)

# ------------------------------------------------------------------------------
# 4. External Dependencies: gbwt & gbwtgraph
# ------------------------------------------------------------------------------
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

set(GBWT_LIBRARY
        "${GBWT_DIR}/lib/libgbwt.a")

set(GBWTGRAPH_LIBRARY
        "${GBWTGRAPH_DIR}/lib/libgbwtgraph.a")

set(GBWT_INCLUDE_DIR
        "${GBWT_DIR}/include")

set(GBWTGRAPH_INCLUDE_DIR
        "${GBWTGRAPH_DIR}/include")

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
set_target_properties(gbwt_lib PROPERTIES
        IMPORTED_LOCATION "${GBWT_LIBRARY}")
add_dependencies(gbwt_lib gbwt_ext)

add_library(gbwtgraph_lib STATIC IMPORTED)
set_target_properties(gbwtgraph_lib PROPERTIES
        IMPORTED_LOCATION "${GBWTGRAPH_LIBRARY}")
add_dependencies(gbwtgraph_lib gbwtgraph_ext)

# ------------------------------------------------------------------------------
# 5. Main Executable Target
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
        src/cdx_writer.cpp
        src/cdx_writer.h
        src/cli.cpp
        src/cli.hpp
)

set_target_properties(cdx_builder PROPERTIES
        BUILD_RPATH "${LIBBDSG_DIR}/lib"
)

target_include_directories(cdx_builder PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${GBWT_INCLUDE_DIR}
        ${GBWTGRAPH_INCLUDE_DIR}
        ${ZSTD_INCLUDE_DIR}
        ${JANSSON_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
)

target_link_libraries(cdx_builder PRIVATE
        libbdsg
        gbwtgraph_lib
        gbwt_lib
        OpenSSL::Crypto
        ${ZSTD_LIBRARY}
        ${JANSSON_LIBRARIES}
        OpenMP::OpenMP_CXX
        Threads::Threads
        CLI11::CLI11
        z
)

# ------------------------------------------------------------------------------
# 6. Google Test
# ------------------------------------------------------------------------------
enable_testing()

add_executable(run_tests
        tests/test_gbz_io.cpp
        src/gbz_IO.cpp
)

set_target_properties(run_tests PROPERTIES
        BUILD_RPATH "${LIBBDSG_DIR}/lib"
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
