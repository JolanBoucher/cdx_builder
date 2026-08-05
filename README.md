# cdx_builder

cdx_builder generates a compact CDX coordinate index from a pangenome graph in GBZ format.

By combining graph topology with haplotype path information, it assigns continuous local
coordinates to nodes within connected components. The resulting CDX index is designed to
enable high-throughput read coverage calculations from GAM alignment files.

---
## Usage

```bash
cdx_builder <input.gbz> [OPTIONS]
```

### Linearization

```text
-i, --iteration INT
    Maximum relaxation iterations.
    Default: 100

-t, --threshold FLOAT
    Convergence threshold.
    Smaller values increase accuracy but may require more iterations.
    Default: 0.01

-l, --lambda-anchor FLOAT
    Balance between path-derived coordinates and topology smoothing.
    Range: [0.0, 1.0]
    Default: 0.7
```

### Output

```text
-o, --output FILE
    Output CDX file path.
    Default: <input>.cdx

-c, --compress LEVEL
    Write a compressed .cdx.zst file.
    Compression level: 1-22
    Default: 3

-d, --debug
    Output a TSV representation to stdout instead of generating a CDX file.
```

### Examples

```bash
# Build a standard CDX index
cdx_builder graph.gbz

# Specify output file
cdx_builder graph.gbz -o graph.cdx

# Compressed output
cdx_builder graph.gbz -c 9

# TSV debug output
cdx_builder graph.gbz -d > debug.tsv

# More accurate relaxation
cdx_builder graph.gbz -t 0.001

# Highly tangled graph with sparse haplotype support
cdx_builder graph.gbz -l 0.3 -i 1000

# Well-supported graph with consistent haplotype structure
cdx_builder graph.gbz -l 0.95

# Display help
cdx_builder -h
```

---

## Features
- **GBZ Native:** Direct loading and validation of GBWT/GBZ pangenome graphs.
- **Haplotype-Aware Indexing:** Incorporates path information to compute stable node coordinates.
- **Component-Level Mapping:** Indexes nodes and metadata per connected component.
- **Flexible Export:** Outputs binary .cdx, Zstandard-compressed .cdx.zst, or TSV format for debugging.

---

## Pipeline at a Glance
1. **Graph Loading:** Reads GBZ graph and extracts node metadata.
2. **Component Analysis:** Partitions graph topology into connected components.
3. **Topology Relaxation:** Optimizes node spatial positioning using path weights and sparse graph layout algorithms.
4. **Coordinate Mapping:** Assigns continuous local coordinates across each component.
5. **Serialization:** Exports the final CDX index.

---

## Dependencies

### Bundled (built automatically by CMake — no manual steps required)

| Dependency    | Purpose                                   |
|---------------|--------------------------------------------|
| SDSL          | Succinct data structures                   |
| GBWT          | Haplotype-aware graph index                |
| GBWTGraph     | Graph representation on top of GBWT        |
| libhandlegraph| Graph interface used by libbdsg            |
| libbdsg       | Pangenome graph backend (GBZ support)      |
| CLI11         | Command-line argument parsing              |
| GoogleTest    | Unit test framework (only when building tests) |
| **cdx_lib**   | CDX format/IO code shared with `cdx_coverage` |

GoogleTest is fetched and built from source automatically via CMake's `FetchContent`
(see `cmake/FetchGoogleTest.cmake`) — no submodule or system package needed.

All of the above (except cdx_lib and GoogleTest) live under `deps/` as git submodules of
this repository.
`cdx_lib` is a **separate repository**, vendored the same way as `deps/cdx_lib` — see
[Getting the code](#getting-the-code) below.

### System dependencies (must be installed separately)

- CMake ≥ 3.16
- Ninja (or Make)
- A C++17-compatible compiler (GCC ≥ 11 or Apple Clang)
- OpenMP runtime
- OpenSSL development libraries
- zstd development libraries
- zlib development libraries
- pkg-config
- Boost
- Jansson
- Git

---

## Getting the code

This repository uses git submodules for every bundled dependency, including `cdx_lib`.
**Clone recursively** to pull everything in one shot:

```bash
git clone --recurse-submodules https://github.com/JolanBoucher/cdx_builder.git
cd cdx_builder
```

If you already have a clone without submodules (or pulled new commits that bump a
submodule), fetch/update them with:

```bash
git submodule update --init --recursive
```

> `cdx_lib` is expected at `deps/cdx_lib`. If you need to build against a local working
> copy instead (e.g. while developing `cdx_lib` itself), point CMake at it directly:
> `cmake -S . -B build -DCDX_LIB_DIR=/path/to/cdx_lib`.

---

## Linux (Ubuntu 20.04+)

### One-command setup

`install_ubuntu20.sh` installs every system package listed below (including a newer
GCC if needed), fetches git submodules, configures and builds the project with
CMake+Ninja, and runs the unit test suite:

```bash
./install_ubuntu20.sh
```

Useful flags: `--no-tests` (build only, skip ctest), `--build-dir <dir>`,
`--build-type <type>`, `--jobs <N>`, `--no-submodules`. Run `./install_ubuntu20.sh --help`
for the full list. The script is safe to re-run.

### Prerequisites (manual setup)

If you'd rather install things yourself instead of using the script above:

```bash
sudo apt update

sudo apt install -y \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    libboost-all-dev \
    libjansson-dev \
    libssl-dev \
    libzstd-dev \
    zlib1g-dev \
    git
```

On GCC, OpenMP support (`libgomp`) ships with `build-essential`, so no extra package is
needed. Ubuntu 20.04's default `cmake` (3.16.3) and `g++` (9.x) both satisfy the minimum
requirements once `build-essential` is up to date; if your GCC is older than 11, install
a newer one (e.g. `sudo apt install g++-11` and select it via `CXX=g++-11`).

### Build

```bash
cmake -S . \
      -B build-linux \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build-linux -j$(nproc)
```

---

## macOS (Apple Silicon / ARM64)

### Prerequisites

Install the required packages using Homebrew:

```bash
brew install \
    libomp \
    openssl@3 \
    zstd \
    cli11 \
    pkg-config \
    boost \
    jansson \
    git
```

Apple Clang does not ship an OpenMP runtime, so `libomp` from Homebrew is required. The
build automatically locates your Homebrew prefix via `brew --prefix` — no manual path
configuration is needed even if Homebrew isn't installed at the default `/opt/homebrew`.

### Build

```bash
cmake -S . \
      -B build-macos \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build-macos
```

---

## Troubleshooting

- **`cdx_lib not found at deps/cdx_lib`** — you cloned without submodules. Run
  `git submodule update --init --recursive` from the repository root.
- **`CMake ... or higher is required` coming from a dependency** — a submodule's
  `CMakeLists.txt` requires a newer CMake than you have installed. Update CMake, or
  check that the submodule is pinned to the expected commit (`git submodule status`).
- **OpenMP not found on macOS** — confirm `brew install libomp` succeeded and that
  `brew --prefix libomp` resolves to a valid path.
- **Stale build after updating a submodule** — remove the build directory
  (`rm -rf build-linux` / `build-macos`) and reconfigure; `ExternalProject`-based
  dependencies (SDSL, GBWT, GBWTGraph) don't always pick up submodule updates in place.

---

## Notes

- All bundled graph libraries (SDSL, GBWT, GBWTGraph, libhandlegraph, libbdsg) are
  compiled automatically during the build process — no dependency-specific `install.sh`
  scripts need to be run manually.
- `cdx_lib` is the only bundled dependency that is *not* built by an `ExternalProject`
  step; it's added directly via `add_subdirectory` from `deps/cdx_lib`.
- OpenMP support is enabled automatically when available.
- The root `CMakeLists.txt` is intentionally thin: compiler checks, macOS/Homebrew
  handling, third-party dependency resolution, and the GoogleTest fetch each live in
  their own file under `cmake/`, included from the root file in dependency order.
- `cmake/FetchGoogleTest.cmake` only fetches GoogleTest if a `gtest` target doesn't
  already exist: libbdsg's vendored copy of sdsl-lite unconditionally builds its own
  bundled googletest snapshot, so on most configurations that vendored copy is reused
  instead of building a second one.
