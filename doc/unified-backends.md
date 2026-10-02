# Native OpenMP and optional RAJA GPU builds

See the [Perlmutter validation record](unified-backends-validation.md) for the
CPU/CUDA regression results and measured four-node performance.

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
CUDA defaults its host compiler to `CMAKE_CXX_COMPILER` before enabling CUDA;
an explicit `CMAKE_CUDA_HOST_COMPILER` or `CUDAHOSTCXX` takes precedence.
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
- `USE_SZ=ON`: SZ plus its HDF5 filter (`H5Z_SZ.h` and filter library),
  located using `SZ_ROOT` and `SZ_FILTER_ROOT`; requires HDF5. Configuration
  checks that the selected filter compiles and links.
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

## Build presets and HIP handoff

CMake 3.21 or newer provides `cpu`, `cuda`, and `hip` presets with independent
build directories. For example:

```sh
cmake --preset cpu -DSW4_BUILD_MOPT=ON
cmake --build --preset cpu -j 16
cmake --preset cuda -DCMAKE_CUDA_ARCHITECTURES=80 \
  -DCMAKE_PREFIX_PATH='/path/to/cuda-raja;/path/to/cuda-umpire'
cmake --build --preset cuda -j 16
cmake --preset hip -DCMAKE_HIP_ARCHITECTURES=gfx90a \
  -DCMAKE_HIP_COMPILER=/path/to/rocm/llvm/bin/clang++ \
  -DCMAKE_CXX_COMPILER=/path/to/rocm/llvm/bin/clang++ \
  -DCMAKE_PREFIX_PATH='/path/to/hip-raja;/path/to/hip-umpire'
cmake --build --preset hip -j 16
```

Load the site's MPI, Fortran, BLAS/LAPACK and optional shared-library modules
first. Select a host compiler and OpenMP runtime compatible with the ROCm
compiler. Host OpenMP defaults to `OFF` for CUDA and `ON` for HIP, matching
their Make workflows. `SW4_GPU_HOST_OPENMP` overrides this without changing
native CPU or inversion targets. For CUDA with host OpenMP enabled, use
`OMP_PROC_BIND=false`: binding a single OpenMP thread to one hardware thread
can starve CUDA progress and make Cartesian refinement appear stalled.
CMake obtains host OpenMP
flags from the compiler instead of assuming GCC flags. `cmake --install`
installs the executables under `CMAKE_INSTALL_PREFIX`.

The HIP setup follows `raja:Makefile.hipcc` and the Frontier configuration:
it compiles the same GPU source manifest, links with the HIP compiler and
uses `-fgpu-rdc` for compilation and linking, with
the platform default stream used by SW4's explicit copies and synchronization.
Build CAMP/RAJA with `-DCAMP_USE_PLATFORM_DEFAULT_STREAM=ON`; the HIP policy
checks the installed configuration, as do the CUDA policies. Modern CAMP rejects a command-line macro
override because it can violate the one-definition rule. `SW4_GPU_MPI_BUFFERS` selects
`STAGED` (host-staged, default for both GPU presets) or `MANAGED` (the legacy
HIP Make setting). Staged validation uses `MPICH_GPU_SUPPORT_ENABLED=0`.
Managed buffers require a compatible memory/MPI transport configuration; the
CUDA curvilinear interface explicitly requires GPU-aware MPI in that mode.
For Cray GPU-aware MPI, set `SW4_GPU_MPI_TRANSPORT_LIBRARY` to the full path of
`libmpi_gtl_cuda.so` or `libmpi_gtl_hsa.so`, as appropriate. CMake retains this
indirectly loaded library in the link and adds its directory to the build RPATH.
It is optional for staged runs with `MPICH_GPU_SUPPORT_ENABLED=0`; native builds
do not discover or link GPU transport libraries. This reproduces the Frontier
Make configuration's HSA transport link without hardcoding a site installation.
HIP builds omit CUDA/NVML and disable legacy roctracer instrumentation. The old
Make config's optional SCR, Caliper and HPCToolkit instrumentation is disabled
by these presets. HIP compilation and runtime acceptance must be
completed on a ROCm system; Perlmutter has no HIP toolchain.

Per-array CUDA prefetch is disabled by default, matching the working RAJA
Make build. Enabling `SW4_CUDA_ARRAY_PREFETCH` is experimental: the legacy
raw-pointer helper targets device zero, and prefetch is outside the accepted
configuration.
`SW4_CUDA_POOL_PREFETCH=ON` separately enables the existing Umpire pool prefetch
before time stepping and requires Umpire. Neither option changes the numerical
kernels. CMake preserves ordinary floating-point compilation; it does not
implicitly enable CUDA fast math.

## Four-node performance acceptance on demand

Export the original profiling input without changing simulation parameters:

```sh
git show raja:performance/large/hmr3.in > "$SCRATCH/hmr3-raja.in"
python tests/backends/profile_hmr3.py \
  --baseline /path/to/exact-raja/sw4 --candidate build/cuda/bin/sw4 \
  --input "$SCRATCH/hmr3-raja.in" \
  --model /path/to/USGSBayAreaVM-08.3.0-corder.rfile \
  --work-dir "$SCRATCH/hmr3-comparison"
```

Run inside an existing allocation of four GPU nodes. The driver launches
16 MPI ranks, four per node, one GPU per rank and one OpenMP thread. It runs
the full nine-second simulation three times per executable in alternating
order. Only the model and output paths are relocated. Logs, executable/input
SHA256 identities, model size/mtime, commands, timings and finite station
waveform comparisons are retained. Acceptance requires the peak-scaled waveform
tolerance used above and a candidate median solver time no more than 5% slower
than the baseline. `--max-slowdown` changes this threshold explicitly. Record
both build configurations alongside the report, especially fast-math and
prefetch differences. `--only baseline --repeats 1` prepares the first run;
`--resume` reuses completed runs only when their identities match.

## Production review follow-up

See [the review correction record](unified-backends-review-fixes.md) for approved
scientific changes, on-demand regression commands and current acceptance limits.
GPU anisotropic operations use the reference host implementation with explicit
stream waits before reading GPU-updated arrays. Anisotropic mesh refinement
remains unsupported by the existing input contract.

GPU event commands share the native name and default-path behavior.
`event parallel=yes` is rejected in GPU builds because existing GPU operations
use world communicators and cannot safely isolate event groups. Native OpenMP
inversion retains parallel-event support. GPU events can run in independent jobs.
This explicit restriction does not establish full GPU inversion/event parity.
