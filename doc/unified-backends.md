# Native OpenMP and optional RAJA GPU builds

The default backend is native OpenMP. It does not discover, include, build, or
link RAJA, Umpire, CUDA, HIP, or NVML. GPU execution is an explicit build choice.
Use separate build directories for different backends and dependency stacks.

```sh
cmake -S . -B build/openmp
cmake --build build/openmp -j

cmake -S . -B build/cuda -DSW4_BACKEND=CUDA \
  -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DRAJA_DIR=/path/to/raja/lib/cmake/raja \
  -Dumpire_DIR=/path/to/umpire/lib64/cmake/umpire
cmake --build build/cuda -j

cmake -S . -B build/hip -DSW4_BACKEND=HIP \
  -DCMAKE_HIP_ARCHITECTURES=gfx90a \
  -DRAJA_DIR=/path/to/hip-raja/lib/cmake/raja \
  -Dumpire_DIR=/path/to/hip-umpire/lib64/cmake/umpire
cmake --build build/hip -j
```

Each configuration produces `bin/sw4`. CUDA/HIP need matching GPU-enabled
RAJA/Umpire packages and compiler toolchains. GPU builds default to C++20 for
the current installed packages; `SW4_GPU_CXX_STANDARD` can select the standard
required by an older dependency installation. CUDA architecture selection
uses CMake's standard `CMAKE_CUDA_ARCHITECTURES`; HIP uses
`CMAKE_HIP_ARCHITECTURES`. HIP configuration requires CMake 3.21 or newer.
The retained RAJA implementation uses double precision. Native builds can
select `SW4_PRECISION=single` as well as the default `double`.

Native material inversion is available with `SW4_BUILD_MOPT=ON`. It produces
`bin/sw4mopt` using native OpenMP, including when the forward `sw4` executable
uses a GPU backend. This does not introduce a GPU inversion implementation.
The existing CPU Makefile remains a native CPU build path; use CMake for the
unified GPU builds.

## Optional shared libraries

All backends use the same switches:

- `USE_PROJ=ON`: PROJ 6 or newer, located using its CMake package.
- `USE_HDF5=ON`: parallel HDF5, located using `HDF5_ROOT` or package defaults.
- `USE_ZFP=ON`: ZFP and H5Z-ZFP CMake packages; requires HDF5.
- `USE_SZ=ON`: SZ headers/library, located using `SZ_ROOT`; requires HDF5.
- `USE_FFTW3=ON`: precision-matched FFTW and MPI FFTW, using `FFTW_ROOT`.
- `SW4_USE_UMPIRE=OFF`: GPU builds use the retained non-Umpire allocator.

Set compilers and MPI compiler wrappers explicitly when the site's environment
contains multiple toolchains. On Perlmutter, load the matching parallel HDF5
module; its wrapper may also need `HDF5_CC` and `HDF5_CLINKER` set to the matching
MPI C compiler. Run substantial builds and simulations on a Slurm allocation.

## Shared functionality and numerical implementations

The integration retains the native sources from `fix-sfile-srf` at
`5e6746ce` and the RAJA sources from `raja` at `52b51c92`. Common executable
code is shared; `SW4_USE_RAJA` selects sections that differ between the
implementations. CMake applies this definition only to GPU targets. It is not
a runtime option and is not a request to run RAJA OpenMP on CPU.

The CPU implementation is the reference for shared SW4 behavior unless evidence
shows it is incorrect. Changes to numerical behavior require scientific review.

Native and RAJA source manifests are explicit in `cmake/SW4Sources.cmake`.
Native material inversion sources are retained. GPU-only kernels, policies,
allocation and profiling helpers, and legacy Fortran kernels are retained.
Sfile and receiver HDF5 readers use one shared implementation, including the
approved interface-index clamping and `USEZVALUE` behavior. SSI dataset
lifetime/progress handling is also shared. CPU and GPU numerical kernels retain
their existing implementations; floating-point results need tolerance-based
comparison rather than byte equality.

## Tests on demand

The existing `pytest/` files and case selection are unchanged. Existing CTest
cases remain available when `BUILD_TESTING=ON`; a result check now depends on
its corresponding run. Additional backend checks are not registered with
pytest or CTest and are invoked explicitly.

Use a shared scratch directory for temporary inputs: login-node `/tmp` is not
visible on compute nodes. With an existing Slurm allocation:

```sh
python tests/backends/run_hdf5.py --backend OPENMP \
  --sw4 build/openmp/bin/sw4 --work-dir "$SCRATCH/sw4-tests/openmp"
python tests/backends/run_hdf5.py --backend CUDA \
  --sw4 build/cuda/bin/sw4 --work-dir "$SCRATCH/sw4-tests/cuda"

python tests/backends/compare_waveforms.py \
  --cpu build/openmp/bin/sw4 --gpu build/cuda/bin/sw4 \
  --work-dir "$SCRATCH/sw4-tests/comparison"

cmake --build build/cuda --target sw4_cuda_compat
srun -n 1 --gpus-per-task=1 build/cuda/bin/sw4_cuda_compat
```

The HDF5 driver reuses the existing receiver metadata, corrupt-restart,
restart-output-offset and compressed SSI regression scripts. It preserves
logs, rejects skips, and supplies GPU launch bindings without editing those
scripts. Run the existing solver suite separately on each backend. Preserve
inputs, compiler options, and output identity when comparing CPU and GPU runs.

The waveform driver checks inline and Sfile materials, topography station
placement, attenuation, mesh refinement, and text/HDF5 ruptures. It generates
its inputs and retains logs and a JSON comparison report. Select individual
cases with repeated `--case` arguments. Waveforms must be finite and nonzero;
the maximum difference for each component must not exceed
`atol + rtol * max(abs(cpu_waveform))`. Defaults are `rtol=2e-5` and
`atol=1e-12`; these are regression thresholds, not a general scientific
acceptance criterion. Requested and sampled station positions are also checked.
SRF coordinates in the fixture are exactly representable in single precision,
so the HDF5 converter's float32 storage does not move the test sources.
