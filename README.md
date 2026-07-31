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
- **GBZ Native :** Direct loading and validation of GBWT/GBZ pangenome graphs.

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
## Build

The following bundled dependencies are built automatically by CMake:

- libbdsg
- libhandlegraph
- GBWT
- GBWTGraph
- SDSL

No manual compilation of bundled dependencies is required.

The following system dependencies must be installed separately:

- CMake ≥ 3.16
- Ninja
- OpenSSL development libraries
- zstd development libraries
- pkg-config
- Boost
- Jansson
- Git
- A C++17-compatible compiler

---

## Linux (Ubuntu 20.04+)

### Prerequisites

Install the required system packages:

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
    git
```

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

### Build

```bash
cmake -S . \
      -B build-macos \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release

cmake --build build-macos
```

---

## Notes

- All bundled graph libraries are compiled automatically during the build process.
- No dependency-specific `install.sh` scripts need to be executed manually.
- The build system automatically compiles:
    - SDSL
    - GBWT
    - GBWTGraph
    - libhandlegraph
    - libbdsg
- OpenMP support is enabled automatically when available.
