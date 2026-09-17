# Building SW4 for Perlmutter A100 GPUs

This guide records the working Perlmutter GPU build as of 2026-09-17. It
targets NVIDIA A100 GPUs (`sm_80`) and uses:

- CUDA 13.2 from the NERSC module environment;
- Cray parallel HDF5 1.14.3.9;
- RAJA 2026.07.0 and Umpire 2026.07.1 from `third_party/`;
- PROJ 9.9.0 from `third_party/`;
- ZFP 1.0.1 built with 8-bit bitstream words;
- H5Z-ZFP 1.1.1 built against that ZFP and parallel HDF5.

Before upgrading, check the official stable releases for
[RAJA](https://github.com/LLNL/RAJA/releases),
[Umpire](https://github.com/LLNL/Umpire/releases),
[PROJ](https://proj.org/en/stable/download.html),
[ZFP](https://github.com/LLNL/zfp/releases), and
[H5Z-ZFP](https://github.com/LLNL/H5Z-ZFP/releases). Do not assume the version
numbers in this document remain current.

## 1. Load the Perlmutter environment

Run all commands from the SW4 repository root:

```bash
module load cudatoolkit/13.2
module load cray-hdf5-parallel/1.14.3.9
module load cray-fftw

export SW4_SOURCE=$PWD
export CUDA_ROOT=$CUDA_HOME
export HDF5_ROOT=/opt/cray/pe/hdf5-parallel/1.14.3.9/gnu/12.3
export HDF5_DIR=$HDF5_ROOT
export MPICH_ROOT=/opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3
```

Confirm the selected toolchain instead of relying only on module names:

```bash
nvcc --version
rg 'HDF5 Version|Parallel HDF5' "$HDF5_ROOT/lib/libhdf5.settings"
```

The expected CUDA compiler is CUDA 13.2, and `libhdf5.settings` must report
`Parallel HDF5: yes`.

On the system used for this build, Lmod printed the non-fatal message
`/opt/cray/pe/lmod/lmod/init/bash: line 100: ERROR:: command not found` while
still applying the requested environment. Verify `CUDA_HOME`, `FFTW_ROOT`, and
the HDF5 path after loading modules if this message persists.

## 2. Build RAJA and Umpire for A100

The source trees must include their Git submodules. Both current releases
require C++20. Their CUDA architecture must match SW4 (`80`, not `90`).

```bash
/usr/bin/cmake -S "$SW4_SOURCE/third_party/RAJA-v2026.07.0" \
  -B "$SW4_SOURCE/third_party/RAJA-v2026.07.0/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SW4_SOURCE/third_party/RAJA-v2026.07.0/install" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_COMPILER="$CUDA_ROOT/bin/nvcc" \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DENABLE_CUDA=ON \
  -DENABLE_OPENMP=OFF \
  -DENABLE_TESTS=OFF \
  -DENABLE_EXAMPLES=OFF
/usr/bin/cmake --build "$SW4_SOURCE/third_party/RAJA-v2026.07.0/build" -j 8
/usr/bin/cmake --install "$SW4_SOURCE/third_party/RAJA-v2026.07.0/build"

/usr/bin/cmake -S "$SW4_SOURCE/third_party/Umpire-v2026.07.1" \
  -B "$SW4_SOURCE/third_party/Umpire-v2026.07.1/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SW4_SOURCE/third_party/Umpire-v2026.07.1/install" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_CUDA_COMPILER="$CUDA_ROOT/bin/nvcc" \
  -DCMAKE_CUDA_HOST_COMPILER=/usr/bin/g++ \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DENABLE_CUDA=ON \
  -DENABLE_OPENMP=OFF \
  -DENABLE_TESTS=OFF \
  -DENABLE_EXAMPLES=OFF
/usr/bin/cmake --build "$SW4_SOURCE/third_party/Umpire-v2026.07.1/build" -j 8
/usr/bin/cmake --install "$SW4_SOURCE/third_party/Umpire-v2026.07.1/build"
```

Inspect both `CMakeCache.txt` files afterward. They should identify CUDA 13.2,
C++20, and `CMAKE_CUDA_ARCHITECTURES=80`.

## 3. Build ZFP with the H5Z-ZFP-compatible word size

This is the most important non-default setting in the dependency stack:

```text
ZFP_BIT_STREAM_WORD_SIZE=8
```

H5Z-ZFP stores the ZFP stream in an HDF5 byte-oriented filter buffer. Build
this ZFP installation with 8-bit bitstream words so its stream representation
matches the H5Z-ZFP integration used by SW4. A default ZFP build uses a
different word size and must not be reused merely because its version matches.

Configure, build, and install ZFP as follows:

```bash
/usr/bin/cmake -S "$SW4_SOURCE/third_party/zfp-1.0.1" \
  -B "$SW4_SOURCE/third_party/zfp-1.0.1/build" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SW4_SOURCE/third_party/zfp-1.0.1/install" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DBUILD_TESTING=ON \
  -DBUILD_UTILITIES=OFF \
  -DZFP_WITH_OPENMP=OFF \
  -DZFP_BIT_STREAM_WORD_SIZE=8
/usr/bin/cmake --build "$SW4_SOURCE/third_party/zfp-1.0.1/build" -j 8
/usr/bin/cmake --install "$SW4_SOURCE/third_party/zfp-1.0.1/build"
```

Always verify the cached value before configuring H5Z-ZFP:

```bash
rg '^ZFP_BIT_STREAM_WORD_SIZE:STRING=8$' \
  third_party/zfp-1.0.1/build/CMakeCache.txt
```

`testviews` passes with this setting. ZFP 1.0.1's `testzfp` regression harness
explicitly requires a 64-bit bitstream word and therefore refuses this 8-bit
configuration. That refusal is not evidence that the H5Z-ZFP configuration is
broken; use the full H5Z-ZFP tests below as the integration test.

## 4. Build H5Z-ZFP against ZFP and parallel HDF5

The HDF5 wrapper needs an explicit MPI backend in this Perlmutter environment:

```bash
export HDF5_CC="$MPICH_ROOT/bin/mpicc"
export HDF5_CLINKER="$MPICH_ROOT/bin/mpicc"

/usr/bin/cmake -S "$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1" \
  -B "$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1/build-cray" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1/install" \
  -DCMAKE_C_COMPILER="$HDF5_ROOT/bin/h5pcc" \
  -DHDF5_PREFER_PARALLEL=ON \
  -DZFP_DIR="$SW4_SOURCE/third_party/zfp-1.0.1/install/lib64/cmake/zfp" \
  -DFORTRAN_INTERFACE=OFF \
  -DBUILD_TESTING=ON
/usr/bin/cmake --build \
  "$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1/build-cray" -j 8
/usr/bin/cmake --install \
  "$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1/build-cray"
```

H5Z-ZFP 1.1.1's test executables needed the math library on this platform. If
the test link reports unresolved `sqrt` or `exp`, ensure these targets in
`third_party/H5Z-ZFP-1.1.1/test/CMakeLists.txt` link `m`:

```cmake
target_link_libraries(test_write_plugin h5z_zfp_shared m)
target_link_libraries(test_write_lib h5z_zfp_static m)
target_link_libraries(test_error h5z_zfp_static m)
```

Run the complete filter test suite:

```bash
export LD_LIBRARY_PATH="$SW4_SOURCE/third_party/zfp-1.0.1/install/lib64:$HDF5_ROOT/lib:${LD_LIBRARY_PATH}"
/usr/bin/ctest \
  --test-dir "$SW4_SOURCE/third_party/H5Z-ZFP-1.1.1/build-cray" \
  --output-on-failure
```

The validated result for this build was 94 of 94 tests passing.

## 5. Build PROJ

PROJ requires a SQLite library and the `sqlite3` command during configuration.
Perlmutter had a usable system SQLite library but no command on `PATH`, so this
build used SQLite 3.53.4 under `third_party/sqlite-autoconf-3530400/install`.
When rebuilding that SQLite source tree, disable ccache because its default
cache under the home filesystem may not be writable from the build context:

```bash
cd "$SW4_SOURCE/third_party/sqlite-autoconf-3530400"
CCACHE_DISABLE=1 CC=/usr/bin/gcc CXX=/usr/bin/g++ \
  ./configure --prefix="$SW4_SOURCE/third_party/sqlite-autoconf-3530400/install"
make -j 8
make install
cd "$SW4_SOURCE"
```

Then build PROJ:

```bash
SQLITE_PREFIX="$SW4_SOURCE/third_party/sqlite-autoconf-3530400/install"

/usr/bin/cmake -S "$SW4_SOURCE/third_party/proj-9.9.0" \
  -B "$SW4_SOURCE/third_party/proj-9.9.0/build-gcc13" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$SW4_SOURCE/third_party/proj-9.9.0/install" \
  -DCMAKE_INSTALL_RPATH="$SQLITE_PREFIX/lib" \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DEXE_SQLITE3="$SQLITE_PREFIX/bin/sqlite3" \
  -DSQLite3_INCLUDE_DIR="$SQLITE_PREFIX/include" \
  -DSQLite3_LIBRARY="$SQLITE_PREFIX/lib/libsqlite3.so" \
  -DBUILD_TESTING=OFF \
  -DBUILD_APPS=OFF \
  -DBUILD_PROJSYNC=OFF
/usr/bin/cmake --build "$SW4_SOURCE/third_party/proj-9.9.0/build-gcc13" -j 8
/usr/bin/cmake --install "$SW4_SOURCE/third_party/proj-9.9.0/build-gcc13"
```

## 6. Configure and build SW4

`configs/make.inc` is intentionally ignored by Git. Preserve it separately
when making a fresh checkout. The working configuration uses these essential
values:

```make
FC = /usr/bin/gfortran
LINKER = /opt/cray/pe/mpich/9.1.0/ofi/gnu/12.3/bin/mpicxx
CXX = $(PREP) nvcc
SW4ROOT = $(CURDIR)

raja_cuda = yes
umpire = yes
fftw = yes
hdf5 = yes
mpi = yes
zfp = yes
proj_6 = yes
openmp = no

RAJA_HOME   = $(CURDIR)/third_party/RAJA-v2026.07.0/install
UMPIRE_HOME = $(CURDIR)/third_party/Umpire-v2026.07.1/install
PROJ_HOME   = $(CURDIR)/third_party/proj-9.9.0/install
HDF5_HOME   = /opt/cray/pe/hdf5-parallel/1.14.3.9/gnu/12.3
ZFP_HOME    = $(CURDIR)/third_party/zfp-1.0.1/install
H5Z_HOME    = $(CURDIR)/third_party/H5Z-ZFP-1.1.1/install

CUDA_ARCH ?= 80
CUDA_ARCH_FLAGS = -arch=sm_$(CUDA_ARCH)
DLINKFLAGS = $(CUDA_ARCH_FLAGS)
```

The complete working configuration is the local `configs/make.inc`. In
addition to the values above, retain the following link requirements:

- C++20 and `$(CUDA_ARCH_FLAGS)` in `EXTRA_CXX_FLAGS`;
- `$CUDA_HOME/lib64/stubs` for link-time `libnvidia-ml` resolution;
- `-lmpi_gtl_cuda` after `-lmpich`, required by Cray MPICH GPU support;
- RPATH entries for HDF5, ZFP, H5Z-ZFP, PROJ, FFTW, MPI, and LibSci;
- `-lsci_gnu`, the current LibSci linker name.

Build SW4:

```bash
make clean
make -j 8 sw4
```

The binary is written to `optimize/sw4`.

## 7. Validate the finished binary

Confirm that every embedded CUDA image targets A100 and that all dynamic
libraries resolve:

```bash
cuobjdump --list-elf optimize/sw4
ldd optimize/sw4 | rg 'not found|libmpi_gtl_cuda|libproj|libhdf5_parallel|libzfp|libcudart'
readelf -d optimize/sw4 | rg 'NEEDED|RUNPATH|RPATH'
```

Expected CUDA image names end in `.sm_80.cubin`. `ldd` must show no
`not found` entries and should resolve `libmpi_gtl_cuda.so`, the parallel HDF5
library, the local PROJ and ZFP libraries, and `libcudart.so.13`.

Do not use a login-node invocation as the GPU runtime test. SW4 initializes
CUDA early and will correctly report `no CUDA-capable device is detected` on a
login node. The repository test driver now detects `NERSC_HOST=perlmutter` and
uses four MPI ranks on one node with:

```text
srun -N 1 -n 4 -c 32 --gpus-per-task=1 \
  --gpu-bind=map_gpu:0,1,2,3 --cpu-bind=cores
```

This follows NERSC's explicit one-rank-per-GPU mapping for a four-GPU node.
`--gpus-per-task=1` also implies per-task GPU binding in current Slurm, so the
explicit map is primarily deterministic and easy to audit; it should not be
treated as a guarantee that it outperforms every other one-GPU-per-rank binding
without an application-specific benchmark.

Run the complete level-0 suite, including all HDF5 cases, in an interactive
GPU allocation:

```bash
salloc --account=m3354_g --constraint=gpu --qos=interactive \
  --time=00:30:00 --nodes=1 --ntasks=4 --cpus-per-task=32 --gpus=4 \
  --job-name=sw4-pytest-a100 \
  /usr/bin/bash -lc \
  'module load cudatoolkit/13.2 && \
   module load cray-hdf5-parallel/1.14.3.9 && \
   module load cray-fftw && \
   export MPICH_GPU_SUPPORT_ENABLED=1 && \
   cd /global/cfs/cdirs/m3354/perl/sw4/pytest && \
   /global/cfs/cdirs/m3354/perl/conda_eqsim/bin/python -u test_sw4.py \
     --level 0 --mpitasks 4 --ompthreads 32 --sw4_exe_dir optimize \
     --pytest_dir /global/cfs/cdirs/m3354/perl/sw4/pytest \
     --usehdf5 4 --verbose'
```

`pytest/test_sw4.py` is SW4's custom test driver rather than a test module for
the Python `pytest` framework. Run it directly as shown above.

The normal Perlmutter build disables application-level prefetching, so also
build and run the focused managed-memory compatibility test inside the GPU
allocation:

```bash
make cuda-compat-test
srun -N 1 -n 1 -c 32 --gpus-per-task=1 --gpu-bind=map_gpu:0 \
  --cpu-bind=cores optimize/test_cuda_compat
```

This exercises device advice, device prefetch, a GPU write, host prefetch, and
host verification through the CUDA-version compatibility adapter.

For backward-compatibility coverage, the same test source was also compiled
successfully with Perlmutter's CUDA 12.9 headers. The production runtime test
above used CUDA 13.2.

## 8. Known source compatibility changes

The current SW4 tree also contains these required compatibility fixes:

- CUDA 13 uses `cudaMemLocation` for `cudaMemPrefetchAsync` and
  `cudaMemAdvise`; `src/CudaCompat.h` provides version-gated adapters so older
  CUDA signatures remain supported.
- The CUDA device-link rule uses `DLINKFLAGS` instead of a hard-coded `sm_90`.
- The Perlmutter configuration derives CUDA compilation and device-link flags
  from one `CUDA_ARCH` value and preserves optional LTO flags.
- H5Z/ZFP variables accept the `H5Z_HOME` and `ZFP_HOME` names used by the
  Perlmutter configuration.
- `TimeSeries` string accessors return const references. A forced CUDA 13.2,
  GCC 13, C++20 rebuild reproduces an nvcc frontend failure in the GCC
  `std::string` copy constructor when these accessors return by value.
