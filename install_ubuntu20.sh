#!/usr/bin/env bash
#
# install_ubuntu20.sh
#
# One-shot setup for cdx_builder on Ubuntu 20.04+: installs system
# dependencies, fetches git submodules, configures and builds the project
# with CMake+Ninja, and (by default) runs the unit test suite via ctest.
#
# Usage:
#   ./install_ubuntu20.sh [options]
#
# Options:
#   --build-dir <dir>   Build directory to configure/build into (default: build-linux)
#   --build-type <type> CMAKE_BUILD_TYPE (default: Release)
#   --jobs <N>           Parallel build jobs (default: nproc)
#   --no-tests           Skip building/running the unit test suite
#   --no-submodules      Skip 'git submodule update --init --recursive'
#   -h, --help            Show this help and exit
#
# Safe to re-run: apt installs, submodule updates, and the CMake
# configure/build steps are all idempotent.

set -euo pipefail

# ------------------------------------------------------------------------------
# Defaults / argument parsing
# ------------------------------------------------------------------------------

BUILD_DIR="build-linux"
BUILD_TYPE="Release"
JOBS="$(nproc)"
RUN_TESTS=1
UPDATE_SUBMODULES=1

usage() {
    sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --build-dir)
            BUILD_DIR="$2"; shift 2 ;;
        --build-type)
            BUILD_TYPE="$2"; shift 2 ;;
        --jobs)
            JOBS="$2"; shift 2 ;;
        --no-tests)
            RUN_TESTS=0; shift ;;
        --no-submodules)
            UPDATE_SUBMODULES=0; shift ;;
        -h|--help)
            usage; exit 0 ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1 ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

log() { printf '\n\033[1;34m==>\033[0m %s\n' "$1"; }

# ------------------------------------------------------------------------------
# 1. System packages
# ------------------------------------------------------------------------------

log "Installing system dependencies (apt)"

sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    git \
    libssl-dev \
    libzstd-dev \
    zlib1g-dev \
    libomp-dev \
    libboost-all-dev \
    libjansson-dev

# Ubuntu 20.04's default GCC is 9.x, but the project requires GCC >= 11
# (enforced in cmake/CompilerRequirements.cmake). Pull a newer compiler from
# the toolchain PPA only if one isn't already available.
GCC_BIN=""
for candidate in g++-13 g++-12 g++-11; do
    if command -v "$candidate" >/dev/null 2>&1; then
        GCC_BIN="$candidate"
        break
    fi
done

if [[ -z "$GCC_BIN" ]]; then
    log "No GCC >= 11 found; installing gcc-11/g++-11 from ppa:ubuntu-toolchain-r/test"
    sudo apt-get install -y software-properties-common
    sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test
    sudo apt-get update
    sudo apt-get install -y gcc-11 g++-11
    GCC_BIN="g++-11"
fi

CXX_COMPILER="$(command -v "$GCC_BIN")"
CC_COMPILER="$(command -v "${GCC_BIN/g++/gcc}")"
log "Using compiler: $CXX_COMPILER"

# ------------------------------------------------------------------------------
# 2. Git submodules (libbdsg, gbwt, gbwtgraph, sdsl-lite, CLI11, cdx_lib)
# ------------------------------------------------------------------------------

if [[ "$UPDATE_SUBMODULES" -eq 1 ]]; then
    log "Fetching git submodules"
    git submodule update --init --recursive
else
    log "Skipping submodule update (--no-submodules)"
fi

# ------------------------------------------------------------------------------
# 3. Configure
# ------------------------------------------------------------------------------

log "Configuring CMake ($BUILD_TYPE, build dir: $BUILD_DIR)"

cmake -S . -B "$BUILD_DIR" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DCMAKE_C_COMPILER="$CC_COMPILER" \
    -DCMAKE_CXX_COMPILER="$CXX_COMPILER" \
    -DCDX_BUILDER_BUILD_TESTS="$([[ $RUN_TESTS -eq 1 ]] && echo ON || echo OFF)"

# ------------------------------------------------------------------------------
# 4. Build
# ------------------------------------------------------------------------------

log "Building (-j$JOBS)"
cmake --build "$BUILD_DIR" -j"$JOBS"

# ------------------------------------------------------------------------------
# 5. Tests
# ------------------------------------------------------------------------------

if [[ "$RUN_TESTS" -eq 1 ]]; then
    log "Running unit tests (ctest)"
    ctest --test-dir "$BUILD_DIR" --output-on-failure -j"$JOBS"
else
    log "Skipping tests (--no-tests)"
fi

log "Done. Binary: $BUILD_DIR/cdx_builder"
