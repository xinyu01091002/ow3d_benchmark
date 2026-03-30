function run_compare_with_ow3d_kinematics_demo()
%RUN_COMPARE_WITH_OW3D_KINEMATICS_DEMO Compare package output with OW3D kinematics probe data.

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(this_dir);
addpath(fullfile(root_dir, 'src'));
package_setup();

cfg = default_config();
cfg.output_dir = fullfile(root_dir, 'output_example', 'ow3d_kinematics_compare_demo');

case_cfg = struct();
case_cfg.data_root = fullfile(fileparts(root_dir), 'uni initial condition', 'ow3d_kinematics_check3');
case_cfg.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
case_cfg.kinematics_file_id = 1;
case_cfg.probe_x_mode = 'auto_centered';
case_cfg.probe_sigma_mode = 'surface';
case_cfg.Alpha = 1.0;
case_cfg.Akp = 0.06;
case_cfg.kd = 1.0;
case_cfg.kp = 0.0279;

ow3d = extract_ow3d_kinematics_probe_case(case_cfg);
input_data = struct();
input_data.t = ow3d.time(:);
input_data.eta_linear = ow3d.eta.first(:);
input_data.h = case_cfg.kd / case_cfg.kp;
input_data.kp = case_cfg.kp;
input_data.g = cfg.g;
input_data.dt = mean(diff(input_data.t));
input_data.fs = 1 / input_data.dt;
input_data.source_file = sprintf('OW3D_KIN:%s', case_cfg.data_root);

result = compute_from_linear_eta_timeseries(input_data, cfg);
export_result_bundle(result, cfg.output_dir, cfg);
metrics = create_ow3d_kinematics_comparison_figure(cfg.output_dir, ow3d, result); %#ok<NASGU>
phi_metrics = create_ow3d_kinematics_phi_decomposition_figure(cfg.output_dir, ow3d, result); %#ok<NASGU>
save(fullfile(cfg.output_dir, 'ow3d_kinematics_probe_reference.mat'), 'ow3d');
fprintf('OW3D kinematics probe comparison demo completed. Output written to:\n%s\n', cfg.output_dir);
fprintf('Selected probe x-index = %d, x = %.6f m (mode: %s)\n', ...
    ow3d.probe_x_index, ow3d.probe_x, ow3d.probe_info.mode);
end
