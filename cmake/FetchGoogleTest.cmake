# ------------------------------------------------------------------------------
# GoogleTest
# ------------------------------------------------------------------------------
#
# Only included when CDX_BUILDER_BUILD_TESTS is ON (see root CMakeLists.txt).
#
# Fetched and built from source instead of relying on a system package: the
# version/ABI of a system-installed gtest is unpredictable across distros, and
# a source build guarantees the plain (non-namespaced) `gtest`/`gtest_main`
# targets that tests/CMakeLists.txt links against. This also means
# install_ubuntu20.sh doesn't need a libgtest-dev package at all -- one less
# thing to install manually.
#
# GUARD: libbdsg vendors its own copy of sdsl-lite, which in turn
# unconditionally add_subdirectory()'s a bundled googletest snapshot
# (deps/libbdsg/bdsg/deps/sdsl-lite/external/googletest) as part of
# ExternalDeps.cmake's add_subdirectory(deps/libbdsg) -- regardless of
# BUILD_TESTING/SDSL_BUILD_TESTS, and regardless of whether *we* want tests at
# all. That already defines `gtest`/`gtest_main` targets (built EXCLUDE_FROM_ALL,
# but the target names exist all the same). If we then also FetchContent our
# own googletest, CMake fails with "add_library cannot create target 'gtest'
# because another target with the same name already exists" (CMP0002).
# So: only fetch our own copy if those targets aren't already provided.
if(TARGET gtest)
    message(STATUS "gtest/gtest_main already provided (vendored via libbdsg's sdsl-lite copy); "
            "reusing it instead of fetching another googletest.")
else()
    include(FetchContent)

    # Match the parent project's runtime library and skip installing gtest itself.
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.15.2
    )
    FetchContent_MakeAvailable(googletest)
endif()
