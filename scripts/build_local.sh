#!/bin/bash
# Build sedfit on taurus without sudo.
#
# All dependencies (compiler + libraries + python) come from the miniforge
# prefix on the taurus userdata disk, so the build and the bulk fits draw
# only from taurus-local storage:
#   /home/taurus/userdata/fmattig/miniforge3
#
# One-time setup (already done):
#   Miniforge3-Linux-x86_64.sh -b -p /home/taurus/userdata/fmattig/miniforge3
#   mamba install -y gxx_linux-64 cmake pkg-config gsl cfitsio libcurl zlib \
#     nlohmann_json numpy pyarrow astropy dustmaps
set -e
cd "$(dirname "$0")/.."
# The login profile exports LD_LIBRARY_PATH pointing at the broken
# gcc-13.2 under /home/norma/hdawson; it must not leak into link or run.
unset LD_LIBRARY_PATH
PREFIX=/home/taurus/userdata/fmattig/miniforge3
PKG_CONFIG_LIBDIR=$PREFIX/lib/pkgconfig:$PREFIX/share/pkgconfig \
  $PREFIX/bin/cmake -B build -S . \
  -DCMAKE_C_COMPILER=$PREFIX/bin/x86_64-conda-linux-gnu-gcc \
  -DCMAKE_CXX_COMPILER=$PREFIX/bin/x86_64-conda-linux-gnu-g++
$PREFIX/bin/cmake --build build -j"$(nproc)"
