# Linear `eta(t)` To Nonlinear Time Series Package

This package is a mostly self-contained MATLAB workflow inside
[OW3D_benchmark](C:/Research/OW3D_benchmark) for converting a unidirectional
single-point linear free-surface time series `eta(t)` into:

- nonlinear free-surface elevation `eta(t)`
- surface horizontal velocity `u_s(t)`
- surface vertical velocity `w_s(t)`

The core linear-eta workflow is now self-contained inside this folder for
normal use. The optional real-OW3D validation scripts still expect
external OW3D data directories.

It exports four groups of results:

- `linear`
- `superharmonic` using VWA
- `subharmonic` using MF12
- `nonlinear`

## What This Version Does

- Takes a `.mat` input containing `t`, `eta_linear`, `h`, `kp`, and optional `g`
- Reconstructs linear surface `u/w` using MF12 first-order kinematics
- Reconstructs second- and third-order superharmonics in time domain using the current VWA coefficient family
- Reconstructs second-order difference-frequency subharmonics using MF12
- Exports a result bundle, CSV table, and summary figures
- Provides OW3D comparison templates, but does not attempt to auto-validate against OW3D

## Important Assumptions

- First version is `unidirectional` only
- Input is a single-point linear time series, not `eta(x,t)` and not a spatial profile `eta(x)`
- The VWA superharmonic path is implemented directly in the time domain
- The MF12 subharmonic path is built by decomposing `eta(t)` into positive-frequency linear components and evaluating MF12 at a fixed spatial point `x=0`
- The exported surface velocities are the package's default physical-surface outputs, but they should still be scientifically checked case by case before publication-grade use

## Input Format

Your `.mat` file should contain:

```matlab
t
eta_linear
h
kp
g   % optional, defaults to 9.81
```

Conventions:

- `t` must be a uniformly spaced column or row vector
- `eta_linear` must have the same length as `t`
- `h` is water depth in meters
- `kp` is peak wavenumber in rad/m
- `g` is gravity in m/s^2

## Quick Start

From the repo root or any location, first initialize the package path:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\startup_package.m');
```

Then run the built-in demo:

```matlab
run_demo_from_sample
```

Run your own input:

```matlab
result = run_from_user_eta_timeseries('C:\path\to\your_input.mat');
```

With a custom output folder:

```matlab
result = run_from_user_eta_timeseries( ...
    'C:\path\to\your_input.mat', ...
    'C:\path\to\desired_output_folder');
```

## Package Layout

- [run_demo_from_sample.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_demo_from_sample.m)
- [run_from_user_eta_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_from_user_eta_timeseries.m)
- [config/default_config.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/config/default_config.m)
- [deps/README.md](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/README.md)
- [src/compute_from_linear_eta_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/compute_from_linear_eta_timeseries.m)
- [examples/generate_sample_input.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/examples/generate_sample_input.m)
- [cases/bingchen/run_bingchen_test_case_batch.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/run_bingchen_test_case_batch.m)
- [cases/bingchen/generate_test_cases_for_bingchen.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/generate_test_cases_for_bingchen.m)
- [cases/bingchen/plot_bingchen_test_case_overview.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/plot_bingchen_test_case_overview.m)
- [validation/real_ow3d/run_compare_with_ow3d_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_demo.m)
- [validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m)
- [validation_templates/compare_with_ow3d_template.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation_templates/compare_with_ow3d_template.m)

## Output Files

A standard run writes:

- `result_bundle.mat`
- `result_timeseries.csv`
- `summary_timeseries.png`
- `summary_spectra_and_closure.png`
- `debug_summary.mat`

The main result structure contains:

```matlab
result.time
result.linear.eta
result.superharmonic.eta
result.subharmonic.eta
result.nonlinear.eta

result.linear.u_surface
result.superharmonic.u_surface
result.subharmonic.u_surface
result.nonlinear.u_surface

result.linear.w_surface
result.superharmonic.w_surface
result.subharmonic.w_surface
result.nonlinear.w_surface

result.meta
```

`result.total` is still kept as a compatibility alias in the MATLAB struct,
but the preferred outward-facing name is now `result.nonlinear`.

## Test Cases For Bingchen Liu

This package now includes one batch generator for Bingchen Liu-requested
wavemaker-style test cases derived from the shared `wavemaker_visual_shared`
signal shape at `x = 0`.

For this batch, the water depth is fixed to `h = 10 m` to match the shared
script, while `T_p`, `A`, and `k_p` are taken directly from the requested
case list.

The current batch uses these three cases:

| Case | `T_p` (s) | `A` (m) | `k_p` (1/m) |
| --- | ---: | ---: | ---: |
| `case_T6_a035_k0047` | 6 | 0.35 | 0.047 |
| `case_T10_a065_k0068` | 10 | 0.65 | 0.068 |
| `case_T14_a100_k0130` | 14 | 1.00 | 0.13 |

Generate the `.mat` inputs only:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\cases\bingchen\generate_test_cases_for_bingchen.m')
```

Generate the inputs and run the full package on all of them:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\cases\bingchen\run_bingchen_test_case_batch.m')
```

Regenerate only the Bingchen overview figures from existing results:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\cases\bingchen\plot_bingchen_test_case_overview.m')
```

This writes the input `.mat` files to:

- [examples/bingchen_test_cases](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/examples/bingchen_test_cases)

and the corresponding package outputs to:

- [output_example/bingchen_test_cases](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/bingchen_test_cases)

The Bingchen output folder also includes:

- `bingchen_input_eta_overview.png`
- `bingchen_nonlinear_overview.png`

## OW3D Validation Template

This package now includes one real `eta` comparison demo based on the
existing extraction workflow from `C:\Research\VWA\VWA time series`, plus a
generic template for later manual extensions.

Run the built-in OW3D comparison demo:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\validation\real_ow3d\run_compare_with_ow3d_demo.m')
```

Run the OW3D kinematics-probe comparison demo:

```matlab
run('C:\Research\OW3D_benchmark\linear_eta_to_nonlinear_timeseries_package\validation\real_ow3d\run_compare_with_ow3d_kinematics_demo.m')
```

This kinematics demo now chooses the probe location automatically by
default, selecting an `x` position where the extracted linear
`\eta^{(1)}(t)` packet is centered as well as possible within the stored
time record. This makes the fixed-point time-series comparison fairer than
using an arbitrary `x` index.

This demo currently compares:

- OW3D `eta^(1)`
- OW3D `eta^(2+) + eta^(3)`
- OW3D `eta^(2-)`
- OW3D total reconstructed `eta`

against the package outputs generated from the extracted OW3D linear
`eta^(1)` time series.

It writes:

- `ow3d_eta_comparison.png`
- `ow3d_eta_comparison_metrics.mat`
- `ow3d_phi_comparison.png`
- `ow3d_phi_comparison_metrics.mat`
- `ow3d_reference_timeseries.mat`

The kinematics-probe demo writes:

- `ow3d_kinematics_probe_comparison.png`
- `ow3d_kinematics_probe_metrics.mat`
- `ow3d_kinematics_probe_reference.mat`

For case-specific extensions, use:

- [validation_templates/compare_with_ow3d_template.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation_templates/compare_with_ow3d_template.m)

That template is still meant to be edited per case. It shows where to load your OW3D reference signals and how to compare them with the package outputs.

## Known Limits

- no directional support yet
- OW3D demo currently validates `eta` only; this package still does not auto-validate `u/w` because the chosen VWA-path data source does not provide ready-to-use OW3D surface velocity time series
- no claim that every quantity is already tuned for every research case
- the time-domain VWA bridge should be treated as a practical first implementation, not a final scientific statement
