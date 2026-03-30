# Vendored Dependencies

This folder contains the small set of external helper functions that were
copied into the package so the core linear-eta workflow can run without
depending on sibling source directories elsewhere in the repo.

## Included VWA Helpers

- [vwa_G_coeff.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/vwa/vwa_G_coeff.m)
- [vwa_mu.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/vwa/vwa_mu.m)
- [vwa_surface_quantity_coeff.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/vwa/vwa_surface_quantity_coeff.m)

These were copied from the repo-level VWA helper folder:

- `test functions for VWA Opensource`

## Included MF12 Helpers

- [mf12_direct_coefficients.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/mf12/mf12_direct_coefficients.m)
- [mf12_second_subharmonic_kinematics.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/deps/mf12/mf12_second_subharmonic_kinematics.m)

These were copied from the repo-level MF12 source folder:

- `irregularWavesMF12/Source`

## Intent

The intent here is not to mirror the full upstream source trees. It is only
to vendor the small subset of functions currently needed by the package's
core workflow.
