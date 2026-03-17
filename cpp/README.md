# C++ Directional OW3D Generator

This directory contains the standalone C++ directional generator for OceanWave3D initial conditions.

## Commands

From the repository root:

```powershell
cmake -S cpp -B cpp/build -G Ninja -DCMAKE_CXX_COMPILER=g++
cmake --build cpp/build
.\cpp\build\ow3d_directional_generator.exe inspect-config cpp\configs\directional_default.json
.\cpp\build\ow3d_directional_generator.exe generate cpp\configs\directional_default.json
```

Simpler one-command runner:

```powershell
.\cpp\run_generator.ps1
```

On Ubuntu/Linux:

```bash
chmod +x cpp/run_generator.sh
./cpp/run_generator.sh
```

This helper now:
- builds with `OpenMP` enabled
- builds with `FFTW` and `OpenMP` enabled
- uses build directory `cpp/build_fftw_omp16`
- runs with `OMP_NUM_THREADS=16` by default

With an explicit config:

```powershell
.\cpp\run_generator.ps1 -Config cpp\configs\directional_default.json
```

Or from `cmd.exe`:

```bat
cpp\run_generator.bat
cpp\run_generator.bat cpp\configs\directional_default.json
```

On Linux with an explicit config:

```bash
./cpp/run_generator.sh cpp/configs/directional_default.json
```

Or with an explicit thread count:

```bash
./cpp/run_generator.sh cpp/configs/directional_default.json generate 8
```

Override the thread count if needed:

```powershell
.\cpp\run_generator.ps1 -Config cpp\configs\directional_default.json -Threads 8
```

```bat
cpp\run_generator.bat cpp\configs\directional_default.json generate 8
```

Verification against MATLAB-style CSV references:

```powershell
.\cpp\build\ow3d_directional_generator.exe verify-matlab cpp\configs\directional_smoke.json <reference_dir>
```

## Current Scope

- full direct MF12 spectral reconstruction for each requested phase
- directional semi-Gaussian spectrum preset
- OW3D export to `OceanWave3D.init` and `OceanWave3D.inp`
- quick-look visualization PNGs plus field CSV dumps

## Dependency Notes

- `FFTW` can be provided either by the system package manager or by placing a local copy under `cpp/third_party/fftw`.
- `OpenMP` is compiler/runtime support, so it needs to come from the Ubuntu toolchain rather than a vendored repo folder.
- The Linux runner prints a short environment check before building so you can see whether FFTW is being found.

## Configs

- `configs/directional_default.json`
  Mirrors the current large directional MATLAB batch settings.
- `configs/directional_smoke.json`
  Small case for fast local validation.
