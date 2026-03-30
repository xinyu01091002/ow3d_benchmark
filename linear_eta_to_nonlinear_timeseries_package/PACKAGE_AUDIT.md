# Package Audit

This note summarizes what is actively used inside
[linear_eta_to_nonlinear_timeseries_package](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package),
what is demo/validation-specific, and what currently looks optional or
lightly used.

## 1. Main End-User Path

These files make up the normal "give me a linear `eta(t)` and return the
package result" workflow.

- [startup_package.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/startup_package.m)
- [run_from_user_eta_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_from_user_eta_timeseries.m)
- [config/default_config.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/config/default_config.m)
- [src/package_setup.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/package_setup.m)
- [src/load_linear_eta_input.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/load_linear_eta_input.m)
- [src/compute_from_linear_eta_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/compute_from_linear_eta_timeseries.m)
- [src/compute_linear_surface_velocity_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/compute_linear_surface_velocity_timeseries.m)
- [src/compute_vwa_superharmonic_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/compute_vwa_superharmonic_timeseries.m)
- [src/compute_mf12_subharmonic_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/compute_mf12_subharmonic_timeseries.m)
- [src/assemble_total_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/assemble_total_timeseries.m)
- [src/export_result_bundle.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/export_result_bundle.m)
- [src/plot_result_summary.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/plot_result_summary.m)
- [src/time_domain_vwa_quantity.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/time_domain_vwa_quantity.m)
- [src/apply_phase_operator_local.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/apply_phase_operator_local.m)
- [src/build_linear_components_from_time_series.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/build_linear_components_from_time_series.m)
- [src/dispersion_wavenumber_from_omega.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/dispersion_wavenumber_from_omega.m)

If you want the package to feel smaller, this is the minimum logical
cluster to think of as the "real package".

## 2. Demo / Example Path

These are helpful, but not required for the main user workflow.

- [run_demo_from_sample.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/run_demo_from_sample.m)
- [examples/generate_sample_input.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/examples/generate_sample_input.m)
- [examples/sample_linear_eta_input.mat](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/examples/sample_linear_eta_input.mat)

## 3. OW3D Validation Path

These files are only needed when doing the OW3D comparisons. They are not
needed for a normal "user input -> package output" run.

- [validation/real_ow3d/run_compare_with_ow3d_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_demo.m)
- [validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation/real_ow3d/run_compare_with_ow3d_kinematics_demo.m)
- [src/build_input_from_ow3d_linear_timeseries.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/build_input_from_ow3d_linear_timeseries.m)
- [src/extract_ow3d_unidirectional_timeseries_case.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/extract_ow3d_unidirectional_timeseries_case.m)
- [src/extract_ow3d_kinematics_probe_case.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/extract_ow3d_kinematics_probe_case.m)
- [src/four_phase_temporal_separation.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/four_phase_temporal_separation.m)
- [src/read_ow3d_surface_bin.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/read_ow3d_surface_bin.m)
- [src/read_ow3d_kinematics_file.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/read_ow3d_kinematics_file.m)
- [src/create_ow3d_eta_comparison_figure.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/create_ow3d_eta_comparison_figure.m)
- [src/create_ow3d_phi_comparison_figure.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/create_ow3d_phi_comparison_figure.m)
- [src/create_ow3d_kinematics_comparison_figure.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/create_ow3d_kinematics_comparison_figure.m)
- [src/create_ow3d_kinematics_phi_decomposition_figure.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/src/create_ow3d_kinematics_phi_decomposition_figure.m)

## 4. Bingchen-Specific Batch Path

These files are useful for the current collaborator handoff, but are not
part of the generic package core.

- [cases/bingchen/run_bingchen_test_case_batch.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/run_bingchen_test_case_batch.m)
- [cases/bingchen/plot_bingchen_test_case_overview.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/plot_bingchen_test_case_overview.m)
- [cases/bingchen/generate_test_cases_for_bingchen.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/cases/bingchen/generate_test_cases_for_bingchen.m)
- [examples/bingchen_test_cases](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/examples/bingchen_test_cases)
- [output_example/bingchen_test_cases](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/bingchen_test_cases)

## 5. Generated Output, Not Source

These directories are useful deliverables, but they are not code and they
make the package tree feel much larger.

- [output_example/demo_run](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/demo_run)
- [output_example/ow3d_compare_demo](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/ow3d_compare_demo)
- [output_example/ow3d_kinematics_compare_demo](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/ow3d_kinematics_compare_demo)
- [output_example/bingchen_test_cases](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/output_example/bingchen_test_cases)

If the package feels crowded, this output area is a major reason.

## 6. Files That Are Currently Lightly Used Or Manual-Only

These are not dead, but they are not part of the automatic core run path.

- [validation_templates/compare_with_ow3d_template.m](C:/Research/OW3D_benchmark/linear_eta_to_nonlinear_timeseries_package/validation_templates/compare_with_ow3d_template.m)
  - referenced in the README only
  - useful as a manual extension template
  - currently not called by any package script

## 7. No Immediate Deletion Candidates In Core Code

Based on the current internal references, there is no obvious core `.m`
file in `src/` that is fully unused. The source tree is crowded more
because it currently mixes:

- core package logic
- demo/example helpers
- OW3D validation utilities
- Bingchen-specific batch utilities
- generated outputs

than because there is a lot of truly dead code.

## 8. Suggested Future Reorganization

Without deleting anything, the clearest future cleanup would be:

1. Keep the generic package entry points at the package root.
2. Move Bingchen-specific scripts into a dedicated subfolder such as
   `cases/bingchen/`.
3. Move OW3D validation scripts into a dedicated subfolder such as
   `validation/real_ow3d/`.
4. Keep `src/` only for reusable code, not case-specific orchestration.
5. Treat `output_example/` as generated artifacts, not as part of the core
   package layout.

That change would reduce visual clutter much more than deleting one or two
files.
