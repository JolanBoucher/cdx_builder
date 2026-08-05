cmake_minimum_required(VERSION 3.16)
project(cdx_builder CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ------------------------------------------------------------------------------
# 1. Build configuration
# ------------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CompilerRequirements)

# ------------------------------------------------------------------------------
# 2. Dependencies
# ------------------------------------------------------------------------------
#
# See cmake/AppleDependencies.cmake and cmake/ExternalDeps.cmake for details.
# Both define variables/targets consumed below and by tests/CMakeLists.txt.

if(APPLE)
    include(AppleDependencies)
endif()

include(ExternalDeps)

# ------------------------------------------------------------------------------
# 3. Main executable target
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
        src/cli.cpp
        src/cli.hpp
        src/cdx_writer.cpp
        src/cdx_writer.h
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
        cdx_lib
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
# 4. Unit tests (Google Test)
# ------------------------------------------------------------------------------
#
# Each testable module gets its own executable, declared in tests/CMakeLists.txt,
# linked only against the dependencies it actually needs.

option(CDX_BUILDER_BUILD_TESTS "Build the cdx_builder unit test suites" ON)

if(CDX_BUILDER_BUILD_TESTS)
    enable_testing()
    include(FetchGoogleTest)
    add_subdirectory(tests)
endif()
