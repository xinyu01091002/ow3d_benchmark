# Self-Contained Status

## Short Answer

This folder is **much closer to self-contained now**, because the core VWA
and MF12 helper files used by the package have been copied into local
dependency folders.

It **is** self-contained for the core package workflow in the current repo
layout, but the real OW3D validation layer still depends on external case
data directories.

## What Works Internally

These paths were checked after the reorganization and still run:

- [run_demo_from_sample.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_demo_from_sample.m)
- [run_from_user_eta_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_from_user_eta_timeseries.m)
- [cases/bingchen/run_bingchen_test_case_batch.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/run_bingchen_test_case_batch.m)
- [validation/real_ow3d/run_compare_with_ow3d_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_demo.m)

So the internal package routing is in good shape.

## Why It Is Not Fully Self-Contained

### 1. Core computation now uses package-local vendored dependencies

[src/package_setup.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/package_setup.m)
now prefers:

- [deps/vwa](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/vwa)
- [deps/mf12](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/mf12)

for the core helper functions used by the package.

The fallback sibling-directory logic is still kept for compatibility, but
the main package computation no longer requires those external source
folders if the local vendored copies are present.

### 2. Real OW3D validation scripts still depend on external data locations

[validation/real_ow3d/run_compare_with_ow3d_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_demo.m)
uses:

- `C:\Research\VWA\VWA time series\unidirectional\timeseriesdata`

[validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m)
uses OW3D case data outside this folder:

- `..\uni initial condition\ow3d_kinematics_check3`

So the validation layer is not portable as-is.

### 3. Generated outputs are still present inside the folder

The folder also includes example outputs in
[output_example](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example),
which is fine for a working project, but it also means the package is not a
minimal clean distribution folder yet.

## Practical Conclusion

### If your question is:
"Can we keep using this folder inside the current repo?"

Answer: **Yes.**

### If your question is:
"Can I zip just this one folder and send it to someone as a standalone tool?"

Answer: **For the core linear-eta package workflow, almost yes. For the
full validation workflow, not yet.**

## What Would Be Needed To Make It Truly Standalone

1. Treat OW3D validation as optional and clearly separate it from the core
   standalone workflow.
2. Optionally ship a lighter distribution without the large generated
   `output_example` tree.

## Current Recommendation

Right now the folder is best understood as:

- a **mostly self-contained core MATLAB package**
- plus **optional validation scripts that still expect external OW3D data**
