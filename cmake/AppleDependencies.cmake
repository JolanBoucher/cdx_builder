# ------------------------------------------------------------------------------
# macOS-specific dependency resolution
# ------------------------------------------------------------------------------
#
# Only included when APPLE is set (see root CMakeLists.txt). Resolves the active
# Homebrew prefix (instead of hardcoding /opt/homebrew) so this also works on
# Intel Macs (/usr/local) or custom Homebrew installs, and configures OpenMP
# support since Apple Clang doesn't ship an OpenMP runtime of its own.

find_program(HOMEBREW_EXECUTABLE brew)
if(HOMEBREW_EXECUTABLE)
    execute_process(
            COMMAND ${HOMEBREW_EXECUTABLE} --prefix
            OUTPUT_VARIABLE HOMEBREW_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()
if(NOT HOMEBREW_PREFIX)
    set(HOMEBREW_PREFIX "/opt/homebrew") # sensible default: Apple Silicon
    message(WARNING "Homebrew not found on PATH; assuming prefix ${HOMEBREW_PREFIX}")
endif()

execute_process(
        COMMAND ${HOMEBREW_EXECUTABLE} --prefix libomp
        OUTPUT_VARIABLE HOMEBREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
)
if(NOT HOMEBREW_LIBOMP_PREFIX)
    set(HOMEBREW_LIBOMP_PREFIX "${HOMEBREW_PREFIX}/opt/libomp")
endif()

# Apple Clang does not ship OpenMP; use Homebrew's libomp explicitly.
set(OpenMP_C_FLAGS
        "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")
set(OpenMP_CXX_FLAGS
        "-Xpreprocessor -fopenmp -I${HOMEBREW_LIBOMP_PREFIX}/include")

set(OpenMP_C_LIB_NAMES "omp")
set(OpenMP_CXX_LIB_NAMES "omp")

set(OpenMP_omp_LIBRARY
        "${HOMEBREW_LIBOMP_PREFIX}/lib/libomp.dylib")

list(APPEND CMAKE_PREFIX_PATH
        "${HOMEBREW_PREFIX}"
        "${HOMEBREW_PREFIX}/opt/openssl@3"
        "${HOMEBREW_PREFIX}/opt/zstd"
        "${HOMEBREW_PREFIX}/opt/cli11"
        "${HOMEBREW_PREFIX}/opt/jansson"
)

# Make sure the linker finds Homebrew-installed libs.
link_directories("${HOMEBREW_PREFIX}/lib")
