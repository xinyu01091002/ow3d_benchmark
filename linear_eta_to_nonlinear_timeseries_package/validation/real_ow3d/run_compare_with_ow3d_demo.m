function run_compare_with_ow3d_demo()
%RUN_COMPARE_WITH_OW3D_DEMO Compare the package eta output with an OW3D case.

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(this_dir);
addpath(fullfile(root_dir, 'src'));
package_setup();

cfg = default_config();
cfg.output_dir = fullfile(root_dir, 'output_example', 'ow3d_compare_demo');

case_cfg = struct();
case_cfg.data_root = 'C:\Research\VWA\VWA time series\unidirectional\timeseriesdata';
case_cfg.Alpha = 8.0;
case_cfg.Akp = 0.02;
case_cfg.kd = 1.0;
case_cfg.kp = 0.0279;
case_cfg.dt = 0.15;
case_cfg.t_steps = 2800:4:4400;
case_cfg.target_idx = 3800;

ow3d = extract_ow3d_unidirectional_timeseries_case(case_cfg);
input_data = build_input_from_ow3d_linear_timeseries(ow3d, cfg.g);
result = compute_from_linear_eta_timeseries(input_data, cfg);
export_result_bundle(result, cfg.output_dir, cfg);
eta_metrics = create_ow3d_eta_comparison_figure(cfg.output_dir, ow3d, result); %#ok<NASGU>
phi_metrics = create_ow3d_phi_comparison_figure(cfg.output_dir, ow3d, result); %#ok<NASGU>

save(fullfile(cfg.output_dir, 'ow3d_reference_timeseries.mat'), 'ow3d');
fprintf('OW3D eta comparison demo completed. Output written to:\n%s\n', cfg.output_dir);
end
