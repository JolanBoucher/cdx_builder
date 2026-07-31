# cdx_builder
build cdx_index from a pangenome gbz format

## Build

All bundled dependencies are built automatically by CMake:

- libbdsg
- libhandlegraph
- GBWT
- GBWTGraph
- pkg-config
- libboost-all-dev
- libjansson-dev
- g++ 11

No manual compilation of dependencies is required.
### Linux (work on 24.04)
#### Prerequisites
Required system packages:

- CMake >= 3.16
- Ninja
- GCC/G++ with C++17 support
- OpenSSL development libraries
- zstd development libraries
- Git

**Example (Ubuntu/Debian):**

```bash
sudo apt update

sudo apt install -y \
  build-essential \
  cmake \
  ninja-build \
  git \
  libssl-dev \
  libzstd-dev
```

#### Build

```bash
cmake -S . -B build-linux -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build build-linux -j$(nproc)
```

### macOS (ARM64)

#### Required packages:

```bash
brew install \
  libomp \
  openssl@3 \
  zstd \
  cli11
```

#### Build:

```bash
cmake -S . -B build-macos -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build build-macos
```

### Notes

The build system automatically compiles all bundled third-party dependencies.
Running dependency-specific `install.sh` scripts manually is not required.