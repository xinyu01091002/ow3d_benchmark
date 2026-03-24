% quick_extract_ow3d_subharmonic_velocity.m
% Strict second-order difference-frequency/subharmonic horizontal velocity
% check against OW3D, with a first-order eta reconstruction sanity check.

clc;
clear;
close all;

addpath(fullfile(pwd, 'irregularWavesMF12', 'Source'));

CFG = struct();
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check3');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.kinematics_file_id = 1;
CFG.time_index = [];
CFG.default_time_index_from_end = 120;
CFG.lambda = 225;
CFG.gravity = 9.81;
CFG.kp_depth = 0.0279;
CFG.sigma_mode = 'surface'; % 'surface' or 'z_target'
CFG.z_target = 0.0;
CFG.z_target_sweep = [0.0, -0.5, -1.0, -2.0, -5.0];
CFG.exclude_surface_layer = false;
CFG.window_source = 'phase0_eta'; % 'phase0_eta' or 'eta11'
CFG.apply_eta11_bandpass = true;
CFG.subharmonic_cutoff_factor = 3.0;
CFG.subharmonic_transition_factor = 0.35;
CFG.linear_fft_rel_tol = 1e-12;
CFG.keep_all_positive_modes = true;
CFG.linear_energy_keep = 0.99999;
CFG.remove_mean_from_subharmonic = false;
CFG.plot_window_lambda = 5.0;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');
CFG.save_mat = true;

four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

phase_data = cell(1, numel(CFG.phi_shifts_deg));
ep_data = cell(1, numel(CFG.phi_shifts_deg));
for idx = 1:numel(CFG.phi_shifts_deg)
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
    kin_path = resolve_kinematics_path_local(case_folder, CFG.kinematics_file_id);
    phase_data{idx} = read_ow3d_kinematics_full_local(kin_path);
    fprintf('Loaded kinematics: %s\n', kin_path);
end

selected_time_index = resolve_time_index_local( ...
    phase_data{1}.n_times_valid, CFG.time_index, CFG.default_time_index_from_end);
x_vec = phase_data{1}.x(:, 1);
sigma_vec = phase_data{1}.sigma(:);
t_selected = phase_data{1}.t(selected_time_index);
kp = 2 * pi / CFG.lambda;
case_kd = extract_kd_from_case_pattern_local(CFG.folder_pattern);
depth_value = case_kd / CFG.kp_depth;

ep_step = resolve_ep_step_from_kinematics_index_local(phase_data{1}, selected_time_index);
for idx = 1:numel(CFG.phi_shifts_deg)
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
    ep_path = fullfile(case_folder, sprintf('EP_%05d.bin', ep_step));
    ep_data{idx} = read_ow3d_ep_snapshot_local(ep_path);
    fprintf('Loaded EP: %s\n', ep_path);
end

eta_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));
etax_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));
u_phases = zeros(numel(CFG.phi_shifts_deg), numel(sigma_vec), numel(x_vec));
w_phases = zeros(numel(CFG.phi_shifts_deg), numel(sigma_vec), numel(x_vec));
phi_phases = zeros(numel(CFG.phi_shifts_deg), numel(sigma_vec), numel(x_vec));
eta_ep_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));
phi_ep_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));

for idx = 1:numel(CFG.phi_shifts_deg)
    eta_phases(idx, :) = squeeze(phase_data{idx}.eta(selected_time_index, :, 1));
    etax_phases(idx, :) = squeeze(phase_data{idx}.etax(selected_time_index, :, 1));
    u_phases(idx, :, :) = squeeze(phase_data{idx}.u(selected_time_index, :, :, 1));
    w_phases(idx, :, :) = squeeze(phase_data{idx}.w(selected_time_index, :, :, 1));
    phi_phases(idx, :, :) = squeeze(phase_data{idx}.phi(selected_time_index, :, :, 1));

    eta_ep_phases(idx, :) = align_ep_to_kinematics_x_local(ep_data{idx}.eta(:), x_vec).';
    phi_ep_phases(idx, :) = align_ep_to_kinematics_x_local(ep_data{idx}.phi(:), x_vec).';
end

eta_harmonics_kin = reconstruct_harmonics_1d_local(eta_phases, four_phase_coef);
etax_harmonics_kin = reconstruct_harmonics_1d_local(etax_phases, four_phase_coef);
u_harmonics = reconstruct_harmonics_xz_local(u_phases, four_phase_coef);
w_harmonics = reconstruct_harmonics_xz_local(w_phases, four_phase_coef);
phi_harmonics = reconstruct_harmonics_xz_local(phi_phases, four_phase_coef);

eta_harmonics_ep = reconstruct_harmonics_1d_local(eta_ep_phases, four_phase_coef);
phi_harmonics_ep = reconstruct_harmonics_1d_local(phi_ep_phases, four_phase_coef);

eta11_kin = eta_harmonics_kin(1, :).';
eta11_ep = eta_harmonics_ep(1, :).';
phi11_ep = phi_harmonics_ep(1, :).';

if CFG.apply_eta11_bandpass
    eta11_kin = frequency_filtering_1d_local(eta11_kin, x_vec, kp, 1);
    eta11_ep = frequency_filtering_1d_local(eta11_ep, x_vec, kp, 1);
    phi11_ep = frequency_filtering_1d_local(phi11_ep, x_vec, kp, 1);
end

sigma_idx = resolve_sigma_index_local(CFG, phase_data, selected_time_index, depth_value, sigma_vec);
sigma_value = sigma_vec(sigma_idx);
z_mean = mean_layer_z_across_phases_local(phase_data, selected_time_index, sigma_value, depth_value);

subharmonic_cutoff = CFG.subharmonic_cutoff_factor * kp;
subharmonic_transition = CFG.subharmonic_transition_factor * kp;
u_fourphase_sub = squeeze(u_harmonics(4, sigma_idx, :));
w_fourphase_sub = squeeze(w_harmonics(4, sigma_idx, :));
phi_fourphase_sub = squeeze(phi_harmonics(4, sigma_idx, :));

u_subharmonic_ow3d = lowpass_component_local(u_fourphase_sub(:), x_vec, subharmonic_cutoff, subharmonic_transition);
w_subharmonic_ow3d = lowpass_component_local(w_fourphase_sub(:), x_vec, subharmonic_cutoff, subharmonic_transition);
phi_subharmonic_ow3d = lowpass_component_local(phi_fourphase_sub(:), x_vec, subharmonic_cutoff, subharmonic_transition);
phix_subharmonic_ow3d = spectral_derivative_x_local(phi_subharmonic_ow3d, x_vec);

u_decomp = decompose_ow3d_u_local( ...
    sigma_value, ...
    squeeze(u_phases(:, sigma_idx, :)), ...
    squeeze(w_phases(:, sigma_idx, :)), ...
    etax_phases, four_phase_coef, x_vec, kp, subharmonic_cutoff, subharmonic_transition);

if CFG.remove_mean_from_subharmonic
    u_subharmonic_ow3d = u_subharmonic_ow3d - mean(u_subharmonic_ow3d);
    w_subharmonic_ow3d = w_subharmonic_ow3d - mean(w_subharmonic_ow3d);
    phi_subharmonic_ow3d = phi_subharmonic_ow3d - mean(phi_subharmonic_ow3d);
    phix_subharmonic_ow3d = phix_subharmonic_ow3d - mean(phix_subharmonic_ow3d);
end

linear_spec = extract_linear_components_local(eta11_kin, x_vec, CFG);
eta11_reconstructed = reconstruct_first_order_eta_local(linear_spec, x_vec, 0.0);
u11_reconstructed = reconstruct_first_order_u_local(linear_spec, x_vec, z_mean, depth_value, CFG.gravity, 0.0);

strict_diff = compute_strict_difference_frequency_u_local( ...
    linear_spec, x_vec, z_mean, depth_value, CFG.gravity, 0.0);
u2minus_theory = strict_diff.u(:);

eta_metrics_kin = compute_metrics_local(eta11_kin, eta11_reconstructed);
eta_metrics_ep = compute_metrics_local(eta11_ep, eta11_reconstructed);
eta_ep_vs_kin_metrics = compute_metrics_local(eta11_kin, eta11_ep);
u1_metrics = compute_metrics_local(squeeze(u_harmonics(1, sigma_idx, :)), u11_reconstructed);
u2minus_metrics = compute_metrics_local(u_subharmonic_ow3d, u2minus_theory);
phix_vs_u_metrics = compute_metrics_local(u_subharmonic_ow3d, phix_subharmonic_ow3d);
u2minus_negative_metrics = compute_metrics_local(u_subharmonic_ow3d, -u2minus_theory);
u2minus_times4_metrics = compute_metrics_local(u_subharmonic_ow3d, 4 * u2minus_theory);
u2minus_negative_times4_metrics = compute_metrics_local(u_subharmonic_ow3d, -4 * u2minus_theory);
z_target_sweep = compute_subharmonic_depth_sweep_local(CFG.z_target_sweep, phase_data, selected_time_index, ...
    depth_value, sigma_vec, u_harmonics, x_vec, linear_spec, CFG.gravity, ...
    subharmonic_cutoff, subharmonic_transition, CFG.remove_mean_from_subharmonic, CFG.exclude_surface_layer);

x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / CFG.lambda;
[x_plot_shifted, eta11_shifted, u_sub_shifted, u_theory_shifted, w_sub_shifted, phi_sub_shifted, phix_shifted, x_center_window] = ...
    recenter_fields_for_plot_local( ...
        x_plot, eta11_kin, eta_phases, u_subharmonic_ow3d, u2minus_theory, ...
        w_subharmonic_ow3d, phi_subharmonic_ow3d, phix_subharmonic_ow3d, CFG.window_source);

eta11_ep_shifted = shift_with_window_local(eta11_ep, eta_phases, eta11_kin, CFG.window_source);
eta11_reconstructed_shifted = shift_with_window_local(eta11_reconstructed, eta_phases, eta11_kin, CFG.window_source);
u11_reconstructed_shifted = shift_with_window_local(u11_reconstructed, eta_phases, eta11_kin, CFG.window_source);
u11_ow3d_shifted = shift_with_window_local(squeeze(u_harmonics(1, sigma_idx, :)), eta_phases, eta11_kin, CFG.window_source);
u_difference_shifted = u_sub_shifted - u_theory_shifted;

x_limits = [x_center_window - 0.5 * CFG.plot_window_lambda, x_center_window + 0.5 * CFG.plot_window_lambda];
x_limits(1) = max(x_limits(1), x_plot_shifted(1));
x_limits(2) = min(x_limits(2), x_plot_shifted(end));
location_label = build_location_label_local(CFG, z_mean);

fig = figure('Color', 'w', 'Position', [120 60 1500 1100]);
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax = nexttile;
plot(ax, x_plot_shifted, eta11_shifted, 'k-', 'LineWidth', 1.8, 'DisplayName', '\eta_{11} from Kinematics'); hold(ax, 'on');
plot(ax, x_plot_shifted, eta11_ep_shifted, '--', 'Color', [0.12 0.45 0.78], 'LineWidth', 1.4, 'DisplayName', '\eta_{11} from EP');
plot(ax, x_plot_shifted, eta11_reconstructed_shifted, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.4, 'DisplayName', 'Reconstructed from linear spectrum');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
ylabel(ax, '\eta^{(1)} (m)');
title(ax, sprintf('First-order eta consistency check (EP aligned by EP(2:4098), t = %.4f s)', t_selected));
legend(ax, 'Location', 'best');

ax = nexttile;
plot(ax, x_plot_shifted, u11_ow3d_shifted, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D u^{(1)}'); hold(ax, 'on');
plot(ax, x_plot_shifted, u11_reconstructed_shifted, '--', 'Color', [0.12 0.45 0.78], 'LineWidth', 1.5, 'DisplayName', 'Reconstructed u^{(1)}');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
ylabel(ax, 'u^{(1)} (m/s)');
title(ax, sprintf('First-order horizontal velocity at %s', location_label));
legend(ax, 'Location', 'best');

ax = nexttile;
plot(ax, x_plot_shifted, u_sub_shifted, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D u_{sub}'); hold(ax, 'on');
plot(ax, x_plot_shifted, u_theory_shifted, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', 'Strict Appendix A u^{(2-)}');
plot(ax, x_plot_shifted, phix_shifted, ':', 'Color', [0.12 0.45 0.78], 'LineWidth', 1.3, 'DisplayName', '\partial_x \phi_{s,sub}');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
ylabel(ax, 'u^{(2-)} (m/s)');
title(ax, sprintf('Subharmonic horizontal velocity at %s', location_label));
legend(ax, 'Location', 'best');

ax = nexttile;
plot(ax, x_plot_shifted, u_difference_shifted, '-', 'Color', [0.50 0.14 0.68], 'LineWidth', 1.7, 'DisplayName', 'OW3D u_{sub} - theory'); hold(ax, 'on');
plot(ax, x_plot_shifted, w_sub_shifted, '--', 'Color', [0.18 0.57 0.30], 'LineWidth', 1.3, 'DisplayName', 'OW3D w_{sub}');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
xlabel(ax, 'x / \lambda');
ylabel(ax, 'Difference / w_{sub}');
title(ax, 'Residual view');
legend(ax, 'Location', 'best');

if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

png_path = fullfile(CFG.output_dir, 'quick_extract_ow3d_subharmonic_velocity.png');
exportgraphics(fig, png_path, 'Resolution', 300);

fig_decomp = figure('Color', 'w', 'Position', [180 120 1500 520]);
ax = axes(fig_decomp);
plot(ax, x_plot_shifted, shift_with_window_local(u_decomp.subharmonic.total, eta_phases, eta11_kin, CFG.window_source), ...
    'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D u_{20} raw'); hold(ax, 'on');
plot(ax, x_plot_shifted, shift_with_window_local(u_decomp.subharmonic.bare, eta_phases, eta11_kin, CFG.window_source), ...
    '--', 'Color', [0.12 0.45 0.78], 'LineWidth', 1.6, 'DisplayName', 'OW3D u_{20} without chain correction');
plot(ax, x_plot_shifted, shift_with_window_local(u_decomp.subharmonic.chain, eta_phases, eta11_kin, CFG.window_source), ...
    '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.6, 'DisplayName', 'OW3D chain correction');
plot(ax, x_plot_shifted, u_theory_shifted, ':', 'Color', [0.45 0.20 0.65], 'LineWidth', 2.0, 'DisplayName', 'Theory u_{20}');
plot(ax, x_plot_shifted, -4 * u_theory_shifted, '-.', 'Color', [0.10 0.55 0.55], 'LineWidth', 1.8, 'DisplayName', '-4 x Theory u_{20}');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
xlabel(ax, 'x / \lambda');
ylabel(ax, 'u_{20} (m/s)');
title(ax, 'u_{20} decomposition: OW3D raw vs chain-corrected split vs theory');
legend(ax, 'Location', 'best');

png_decomp_path = fullfile(CFG.output_dir, 'quick_extract_ow3d_u20_decomposition.png');
exportgraphics(fig_decomp, png_decomp_path, 'Resolution', 300);

fig_scale = figure('Color', 'w', 'Position', [220 150 1500 520]);
ax = axes(fig_scale);
plot(ax, x_plot_shifted, shift_with_window_local(u_decomp.subharmonic.bare, eta_phases, eta11_kin, CFG.window_source), ...
    'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D u_{20} without chain correction'); hold(ax, 'on');
plot(ax, x_plot_shifted, 4 * u_theory_shifted, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', '+4 x Theory u_{20}');
plot(ax, x_plot_shifted, -4 * u_theory_shifted, '-.', 'Color', [0.10 0.55 0.55], 'LineWidth', 1.8, 'DisplayName', '-4 x Theory u_{20}');
xline(ax, x_center_window, '--', 'Color', [0.65 0.65 0.65], 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
grid(ax, 'on');
box(ax, 'on');
xlim(ax, x_limits);
xlabel(ax, 'x / \lambda');
ylabel(ax, 'u_{20} (m/s)');
title(ax, 'OW3D u_{20} without chain correction vs \pm 4 x theory u_{20}');
legend(ax, 'Location', 'best');

png_scale_path = fullfile(CFG.output_dir, 'quick_extract_ow3d_u20_scaling_sanity.png');
exportgraphics(fig_scale, png_scale_path, 'Resolution', 300);

u20_bare_shifted = shift_with_window_local(u_decomp.subharmonic.bare, eta_phases, eta11_kin, CFG.window_source);
[k_abs_plot, amp_bare, amp_theory_p4, amp_theory_m4] = compute_single_sided_spectrum_local( ...
    u20_bare_shifted, 4 * u_theory_shifted, -4 * u_theory_shifted, x_vec);

fig_spec = figure('Color', 'w', 'Position', [240 180 1500 560]);
ax = axes(fig_spec);
plot(ax, k_abs_plot / kp, amp_bare, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D u_{20} without chain correction'); hold(ax, 'on');
plot(ax, k_abs_plot / kp, amp_theory_p4, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.8, 'DisplayName', '+4 x Theory u_{20}');
plot(ax, k_abs_plot / kp, amp_theory_m4, '-.', 'Color', [0.10 0.55 0.55], 'LineWidth', 1.8, 'DisplayName', '-4 x Theory u_{20}');
xline(ax, subharmonic_cutoff / kp, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'DisplayName', 'Subharmonic cutoff');
grid(ax, 'on');
box(ax, 'on');
xlabel(ax, '|k| / k_p');
ylabel(ax, 'Single-sided amplitude');
title(ax, 'Spectrum of OW3D u_{20} without chain correction vs \pm 4 x theory u_{20}');
legend(ax, 'Location', 'northeast');

png_spec_path = fullfile(CFG.output_dir, 'quick_extract_ow3d_u20_scaling_spectrum.png');
exportgraphics(fig_spec, png_spec_path, 'Resolution', 300);

results = struct();
results.x = x_vec(:);
results.x_plot = x_plot(:);
results.x_plot_shifted = x_plot_shifted(:);
results.time_index = selected_time_index;
results.time_value = t_selected;
results.ep_step = ep_step;
results.depth_value = depth_value;
results.kp = kp;
results.sigma_idx = sigma_idx;
results.sigma_value = sigma_value;
results.z_mean = z_mean;
results.linear_spec = linear_spec;
results.strict_difference = strict_diff;
results.eta11_kin = eta11_kin(:);
results.eta11_ep = eta11_ep(:);
results.phi11_ep = phi11_ep(:);
results.eta11_reconstructed = eta11_reconstructed(:);
results.u11_ow3d = squeeze(u_harmonics(1, sigma_idx, :));
results.u11_reconstructed = u11_reconstructed(:);
results.u_subharmonic_ow3d = u_subharmonic_ow3d(:);
results.w_subharmonic_ow3d = w_subharmonic_ow3d(:);
results.phi_subharmonic_ow3d = phi_subharmonic_ow3d(:);
results.phix_subharmonic_ow3d = phix_subharmonic_ow3d(:);
results.ow3d_u_decomposition = u_decomp;
results.u2minus_theory = u2minus_theory(:);
results.metrics = struct();
results.metrics.eta_reconstructed_vs_kin = eta_metrics_kin;
results.metrics.eta_reconstructed_vs_ep = eta_metrics_ep;
results.metrics.eta_ep_vs_kin = eta_ep_vs_kin_metrics;
results.metrics.u1_theory_vs_ow3d = u1_metrics;
results.metrics.u2minus_theory_vs_ow3d = u2minus_metrics;
results.metrics.neg_u2minus_theory_vs_ow3d = u2minus_negative_metrics;
results.metrics.u2minus_times4_theory_vs_ow3d = u2minus_times4_metrics;
results.metrics.neg_u2minus_times4_theory_vs_ow3d = u2minus_negative_times4_metrics;
results.metrics.phixsub_vs_ow3d_usub = phix_vs_u_metrics;
results.metrics.ow3d_u_decomposition = u_decomp.metrics;
results.metrics.depth_sweep = z_target_sweep;
results.window_source = CFG.window_source;
results.x_center_window = x_center_window;
results.cutoff = subharmonic_cutoff;
results.transition = subharmonic_transition;

if CFG.save_mat
    mat_path = fullfile(CFG.output_dir, 'quick_extract_ow3d_subharmonic_velocity.mat');
    save(mat_path, 'results', 'CFG');
else
    mat_path = '';
end

fprintf('\n=== First-order eta reconstruction check ===\n');
fprintf('Kinematics eta11 vs reconstructed: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    eta_metrics_kin.corr, eta_metrics_kin.rmse, eta_metrics_kin.peak_ratio);
fprintf('EP eta11 vs reconstructed:         corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    eta_metrics_ep.corr, eta_metrics_ep.rmse, eta_metrics_ep.peak_ratio);
fprintf('EP eta11 vs Kinematics eta11:      corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    eta_ep_vs_kin_metrics.corr, eta_ep_vs_kin_metrics.rmse, eta_ep_vs_kin_metrics.peak_ratio);

fprintf('\n=== First-order u check at selected z ===\n');
fprintf('OW3D u^(1) vs reconstructed:       corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u1_metrics.corr, u1_metrics.rmse, u1_metrics.peak_ratio);

fprintf('\n=== Strict difference-frequency horizontal velocity ===\n');
fprintf('OW3D u_sub vs strict u^(2-):       corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u2minus_metrics.corr, u2minus_metrics.rmse, u2minus_metrics.peak_ratio);
fprintf('OW3D u_sub vs -strict u^(2-):      corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u2minus_negative_metrics.corr, u2minus_negative_metrics.rmse, u2minus_negative_metrics.peak_ratio);
fprintf('OW3D u_sub vs 4*strict u^(2-):     corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u2minus_times4_metrics.corr, u2minus_times4_metrics.rmse, u2minus_times4_metrics.peak_ratio);
fprintf('OW3D u_sub vs -4*strict u^(2-):    corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u2minus_negative_times4_metrics.corr, u2minus_negative_times4_metrics.rmse, u2minus_negative_times4_metrics.peak_ratio);
fprintf('OW3D u_sub vs d/dx(phi_s,sub):     corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    phix_vs_u_metrics.corr, phix_vs_u_metrics.rmse, phix_vs_u_metrics.peak_ratio);
fprintf('Strict difference-frequency pairs retained: %d\n', strict_diff.num_pairs);

fprintf('\n=== OW3D u decomposition at selected level ===\n');
fprintf('Subharmonic: u_sub vs bare phi_x(sigma):      corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u_decomp.metrics.subharmonic.bare_vs_total.corr, ...
    u_decomp.metrics.subharmonic.bare_vs_total.rmse, ...
    u_decomp.metrics.subharmonic.bare_vs_total.peak_ratio);
fprintf('Subharmonic: u_sub vs chain correction:       corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u_decomp.metrics.subharmonic.chain_vs_total.corr, ...
    u_decomp.metrics.subharmonic.chain_vs_total.rmse, ...
    u_decomp.metrics.subharmonic.chain_vs_total.peak_ratio);
fprintf('Superharmonic: u_2 vs bare phi_x(sigma):      corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u_decomp.metrics.superharmonic.bare_vs_total.corr, ...
    u_decomp.metrics.superharmonic.bare_vs_total.rmse, ...
    u_decomp.metrics.superharmonic.bare_vs_total.peak_ratio);
fprintf('Superharmonic: u_2 vs chain correction:       corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    u_decomp.metrics.superharmonic.chain_vs_total.corr, ...
    u_decomp.metrics.superharmonic.chain_vs_total.rmse, ...
    u_decomp.metrics.superharmonic.chain_vs_total.peak_ratio);

fprintf('\n=== Fixed-z sweep for OW3D u_sub vs strict u^(2-) ===\n');
fprintf(' target_z(m)   sigma_idx    mean_z(m)        corr          rmse        peak_ratio\n');
for ii = 1:numel(z_target_sweep)
    row = z_target_sweep(ii);
    fprintf('%11.3f %10d %12.6f %12.6f %12.6e %12.6f\n', ...
        row.target_z, row.sigma_idx, row.z_mean, row.metrics.corr, row.metrics.rmse, row.metrics.peak_ratio);
end

fprintf('\nSaved figure: %s\n', png_path);
fprintf('Saved decomposition figure: %s\n', png_decomp_path);
fprintf('Saved scaling sanity figure: %s\n', png_scale_path);
fprintf('Saved scaling spectrum figure: %s\n', png_spec_path);
if ~isempty(mat_path)
    fprintf('Saved MAT: %s\n', mat_path);
end
fprintf('Time index = %d, t = %.6f s, EP step = %d\n', selected_time_index, t_selected, ep_step);
fprintf('Depth = %.6f m from kd = %.4f and kp = %.4f 1/m\n', depth_value, case_kd, CFG.kp_depth);
fprintf('Chosen sigma index = %d, sigma = %.6f, mean z = %.6f m\n', sigma_idx, sigma_value, z_mean);
fprintf('Window source = %s, recentered plot peak at x/lambda = %.6f\n', CFG.window_source, x_center_window);
fprintf('Subharmonic cutoff = %.6f 1/m, transition = %.6f 1/m\n', subharmonic_cutoff, subharmonic_transition);
fprintf('max|u_sub_ow3d| = %.6e, max|u2minus_theory| = %.6e, max|w_sub_ow3d| = %.6e\n', ...
    max(abs(u_subharmonic_ow3d)), max(abs(u2minus_theory)), max(abs(w_subharmonic_ow3d)));

function kin_path = resolve_kinematics_path_local(case_folder, file_id)
    if file_id < 10
        file_name = sprintf('Kinematics0%d.bin', file_id);
    else
        file_name = sprintf('Kinematics%d.bin', file_id);
    end
    kin_path = fullfile(case_folder, file_name);
end

function time_index = resolve_time_index_local(n_times_valid, requested_index, default_from_end)
    if isempty(requested_index)
        time_index = max(1, n_times_valid - default_from_end + 1);
        return;
    end

    if requested_index > 0
        time_index = min(requested_index, n_times_valid);
        return;
    end

    time_index = max(1, n_times_valid + requested_index + 1);
end

function data = read_ow3d_kinematics_full_local(file_path)
    fid = fopen(file_path, 'r', 'ieee-le');
    if fid < 0
        error('Could not open kinematics file: %s', file_path);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fread(fid, 1, 'int32');
    xbeg = fread(fid, 1, 'int32');
    xend = fread(fid, 1, 'int32');
    xstride = fread(fid, 1, 'int32');
    ybeg = fread(fid, 1, 'int32');
    yend = fread(fid, 1, 'int32');
    ystride = fread(fid, 1, 'int32');
    tbeg = fread(fid, 1, 'int32');
    tend = fread(fid, 1, 'int32');
    tstride = fread(fid, 1, 'int32');
    dt = fread(fid, 1, 'double');
    nz = fread(fid, 1, 'int32');
    sigma = zeros(nz, 1);
    fread(fid, 2, 'int32');

    nx = floor((xend - xbeg) / xstride) + 1;
    ny = floor((yend - ybeg) / ystride) + 1;
    nt = floor((tend - tbeg) / tstride) + 1;

    tmp = zeros(nx * ny * max(nz, 5), 1);
    tmp(1:5 * nx * ny) = fread(fid, 5 * nx * ny, 'double');
    fread(fid, 2, 'int32');

    x = zeros(nx, ny);
    y = zeros(nx, ny);
    h = zeros(nx, ny);
    x(:) = tmp(1:5:5 * nx * ny);
    y(:) = tmp(2:5:5 * nx * ny);
    h(:) = tmp(3:5:5 * nx * ny);

    for i = 1:nz
        sigma(i) = fread(fid, 1, 'double');
    end
    fread(fid, 2, 'int32');

    eta = zeros(nt, nx, ny);
    etax = zeros(nt, nx, ny);
    phi = zeros(nt, nz, nx, ny);
    u = zeros(nt, nz, nx, ny);
    v = zeros(nt, nz, nx, ny); %#ok<NASGU>
    w = zeros(nt, nz, nx, ny);
    uz = zeros(nt, nz, nx, ny); %#ok<NASGU>
    vz = zeros(nt, nz, nx, ny); %#ok<NASGU>
    wz = zeros(nt, nz, nx, ny); %#ok<NASGU>
    t = (0:nt - 1) * dt * tstride;

    it_valid = 0;
    for it = 1:nt - 1
        tmp_eta = fread(fid, nx * ny, 'double');
        if numel(tmp_eta) < nx * ny
            break;
        end
        eta(it, :) = tmp_eta;
        fread(fid, 2, 'int32');

        tmp_etax = fread(fid, nx * ny, 'double');
        if numel(tmp_etax) < nx * ny
            break;
        end
        etax(it, :) = tmp_etax;
        fread(fid, 2, 'int32');

        tmp_etay = fread(fid, nx * ny, 'double'); %#ok<NASGU>
        if numel(tmp_etay) < nx * ny
            break;
        end
        fread(fid, 2, 'int32');

        tmp_phi = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_phi) < nx * ny * nz
            break;
        end
        phi(it, :) = tmp_phi;
        fread(fid, 2, 'int32');

        tmp_u = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_u) < nx * ny * nz
            break;
        end
        u(it, :) = tmp_u;
        fread(fid, 2, 'int32');

        tmp_v = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_v) < nx * ny * nz
            break;
        end
        v(it, :) = tmp_v;
        fread(fid, 2, 'int32');

        tmp_w = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_w) < nx * ny * nz
            break;
        end
        w(it, :) = tmp_w;
        fread(fid, 2, 'int32');

        tmp_wz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_wz) < nx * ny * nz
            break;
        end
        wz(it, :) = tmp_wz;
        fread(fid, 2, 'int32');

        tmp_uz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_uz) < nx * ny * nz
            break;
        end
        uz(it, :) = tmp_uz;
        fread(fid, 2, 'int32');

        tmp_vz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_vz) < nx * ny * nz
            break;
        end
        vz(it, :) = tmp_vz;
        fread(fid, 2, 'int32');

        it_valid = it;
    end

    if it_valid <= 0
        error('No complete stored kinematics time step could be read from %s', file_path);
    end

    eta = eta(1:it_valid, :, :);
    etax = etax(1:it_valid, :, :);
    phi = phi(1:it_valid, :, :, :);
    u = u(1:it_valid, :, :, :);
    w = w(1:it_valid, :, :, :);
    t = t(1:it_valid);

    data.eta = eta;
    data.etax = etax;
    data.phi = phi;
    data.u = u;
    data.w = w;
    data.x = x;
    data.y = y;
    data.h = h;
    data.sigma = sigma;
    data.t = t;
    data.n_times_valid = it_valid;
    data.dt = dt;
    data.tbeg = tbeg;
    data.tend = tend;
    data.tstride = tstride;
end

function data = read_ow3d_ep_snapshot_local(file_path)
    fid = fopen(file_path, 'r', 'ieee-le');
    if fid < 0
        error('Could not open EP file: %s', file_path);
    end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

    fread(fid, 1, 'int32');
    nx = fread(fid, 1, 'int32');
    ny = fread(fid, 1, 'int32');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    x = fread(fid, [nx ny], 'float64');
    y = fread(fid, [nx ny], 'float64');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    eta = fread(fid, [nx ny], 'float64');
    phi = fread(fid, [nx ny], 'float64');

    data.x = x;
    data.y = y;
    data.eta = eta;
    data.phi = phi;
end

function ep_step = resolve_ep_step_from_kinematics_index_local(kin_data, time_index)
    ep_step = (kin_data.tbeg - 1) + (time_index - 1) * kin_data.tstride;
end

function aligned = align_ep_to_kinematics_x_local(ep_field, x_kin)
    ep_field = ep_field(:);
    x_kin = x_kin(:);

    if numel(ep_field) == numel(x_kin)
        aligned = ep_field;
        return;
    end

    if numel(ep_field) == numel(x_kin) + 2
        aligned = ep_field(2:end - 1);
        return;
    end

    error('Unexpected EP/Kinematics x-size mismatch: EP has %d points, Kinematics has %d points.', ...
        numel(ep_field), numel(x_kin));
end

function harmonics = reconstruct_harmonics_1d_local(fields_by_phase, coef)
    all_fields = [
        fields_by_phase(1, :);
        fields_by_phase(2, :);
        fields_by_phase(3, :);
        fields_by_phase(4, :);
        fields_by_phase(1, :).^2;
        fields_by_phase(2, :).^2;
        fields_by_phase(3, :).^2;
        fields_by_phase(4, :).^2];

    harmonics = zeros(4, size(fields_by_phase, 2));
    for n = 1:4
        harmonics(n, :) = coef(n, :) * all_fields;
    end
end

function harmonics = reconstruct_harmonics_xz_local(fields_by_phase, coef)
    n_phases = size(fields_by_phase, 1);
    nz = size(fields_by_phase, 2);
    nx = size(fields_by_phase, 3);
    if n_phases ~= 4
        error('Expected exactly four phases for four-phase separation.');
    end

    harmonics = zeros(4, nz, nx);
    for iz = 1:nz
        fields_1d = squeeze(fields_by_phase(:, iz, :));
        all_fields = [
            fields_1d(1, :);
            fields_1d(2, :);
            fields_1d(3, :);
            fields_1d(4, :);
            fields_1d(1, :).^2;
            fields_1d(2, :).^2;
            fields_1d(3, :).^2;
            fields_1d(4, :).^2];
        for n = 1:4
            harmonics(n, iz, :) = coef(n, :) * all_fields;
        end
    end
end

function sigma_idx = resolve_sigma_index_local(CFG, phase_data, time_index, depth_value, sigma_vec)
    switch lower(CFG.sigma_mode)
        case 'surface'
            sigma_idx = numel(sigma_vec);
        case 'z_target'
            sigma_idx = choose_sigma_index_near_target_z_local(phase_data, time_index, depth_value, ...
                sigma_vec, CFG.z_target, CFG.exclude_surface_layer);
        otherwise
            error('Unsupported CFG.sigma_mode: %s. Use ''surface'' or ''z_target''.', CFG.sigma_mode);
    end
end

function sigma_idx = choose_sigma_index_near_target_z_local(phase_data, time_index, depth_value, sigma_vec, z_target, exclude_surface)
    z_cost = inf(numel(sigma_vec), 1);
    surface_idx = find(abs(sigma_vec - max(sigma_vec)) < 1e-12, 1, 'first');

    for j = 1:numel(sigma_vec)
        if exclude_surface && j == surface_idx
            continue;
        end

        z_samples = [];
        for p = 1:numel(phase_data)
            eta_now = squeeze(phase_data{p}.eta(time_index, :, 1));
            z_now = -depth_value + sigma_vec(j) .* (depth_value + eta_now(:));
            z_samples = [z_samples; z_now(:)]; %#ok<AGROW>
        end
        z_cost(j) = mean(abs(z_samples - z_target));
    end

    [~, sigma_idx] = min(z_cost);
end

function z_mean = mean_layer_z_across_phases_local(phase_data, time_index, sigma_value, depth_value)
    z_samples = [];
    for p = 1:numel(phase_data)
        eta_now = squeeze(phase_data{p}.eta(time_index, :, 1));
        z_now = -depth_value + sigma_value .* (depth_value + eta_now(:));
        z_samples = [z_samples; z_now(:)]; %#ok<AGROW>
    end
    z_mean = mean(z_samples);
end

function sweep = compute_subharmonic_depth_sweep_local(z_targets, phase_data, time_index, depth_value, sigma_vec, ...
        u_harmonics, x_vec, linear_spec, gravity, subharmonic_cutoff, subharmonic_transition, remove_mean, exclude_surface_layer)
    sweep = repmat(struct( ...
        'target_z', 0, ...
        'sigma_idx', 0, ...
        'sigma_value', 0, ...
        'z_mean', 0, ...
        'u_subharmonic_ow3d', [], ...
        'u2minus_theory', [], ...
        'metrics', struct('corr', NaN, 'rmse', NaN, 'peak_ratio', NaN)), 1, numel(z_targets));

    local_cfg = struct();
    local_cfg.sigma_mode = 'z_target';
    local_cfg.exclude_surface_layer = exclude_surface_layer;

    for ii = 1:numel(z_targets)
        local_cfg.z_target = z_targets(ii);
        sigma_idx = resolve_sigma_index_local(local_cfg, phase_data, time_index, depth_value, sigma_vec);
        sigma_value = sigma_vec(sigma_idx);
        z_mean = mean_layer_z_across_phases_local(phase_data, time_index, sigma_value, depth_value);

        u_ow3d = lowpass_component_local(squeeze(u_harmonics(4, sigma_idx, :)), x_vec, subharmonic_cutoff, subharmonic_transition);
        if remove_mean
            u_ow3d = u_ow3d - mean(u_ow3d);
        end

        strict_diff = compute_strict_difference_frequency_u_local(linear_spec, x_vec, z_mean, depth_value, gravity, 0.0);
        u_theory = strict_diff.u(:);
        metrics = compute_metrics_local(u_ow3d, u_theory);

        sweep(ii).target_z = z_targets(ii);
        sweep(ii).sigma_idx = sigma_idx;
        sweep(ii).sigma_value = sigma_value;
        sweep(ii).z_mean = z_mean;
        sweep(ii).u_subharmonic_ow3d = u_ow3d(:);
        sweep(ii).u2minus_theory = u_theory(:);
        sweep(ii).metrics = metrics;
    end
end

function out = decompose_ow3d_u_local(sigma_value, u_phase_slice, w_phase_slice, etax_phases, coef, x_vec, kp, subharmonic_cutoff, subharmonic_transition)
    u_phase_slice = squeeze(u_phase_slice);
    w_phase_slice = squeeze(w_phase_slice);
    etax_phases = squeeze(etax_phases);

    chain_phase = -sigma_value .* etax_phases .* w_phase_slice;
    bare_phase = u_phase_slice - chain_phase;

    bare_h = reconstruct_harmonics_1d_local(bare_phase, coef);
    chain_h = reconstruct_harmonics_1d_local(chain_phase, coef);
    total_h = reconstruct_harmonics_1d_local(u_phase_slice, coef);

    total_super = frequency_filtering_1d_local(total_h(2, :).', x_vec, kp, 2);
    bare_super = frequency_filtering_1d_local(bare_h(2, :).', x_vec, kp, 2);
    chain_super = frequency_filtering_1d_local(chain_h(2, :).', x_vec, kp, 2);

    total_sub = lowpass_component_local(total_h(4, :).', x_vec, subharmonic_cutoff, subharmonic_transition);
    bare_sub = lowpass_component_local(bare_h(4, :).', x_vec, subharmonic_cutoff, subharmonic_transition);
    chain_sub = lowpass_component_local(chain_h(4, :).', x_vec, subharmonic_cutoff, subharmonic_transition);

    out = struct();
    out.sigma_value = sigma_value;
    out.total_harmonics = total_h;
    out.bare_harmonics = bare_h;
    out.chain_harmonics = chain_h;
    out.superharmonic = struct('total', total_super(:), 'bare', bare_super(:), 'chain', chain_super(:));
    out.subharmonic = struct('total', total_sub(:), 'bare', bare_sub(:), 'chain', chain_sub(:));
    out.metrics = struct();
    out.metrics.subharmonic = struct( ...
        'bare_vs_total', compute_metrics_local(total_sub, bare_sub), ...
        'chain_vs_total', compute_metrics_local(total_sub, chain_sub), ...
        'sum_vs_total', compute_metrics_local(total_sub, bare_sub + chain_sub));
    out.metrics.superharmonic = struct( ...
        'bare_vs_total', compute_metrics_local(total_super, bare_super), ...
        'chain_vs_total', compute_metrics_local(total_super, chain_super), ...
        'sum_vs_total', compute_metrics_local(total_super, bare_super + chain_super));
end

function out = lowpass_component_local(field_in, x_vec, k_cutoff, transition)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = exp(-(abs(kx) / max(transition, dkx)).^4);
    mask(abs(kx) <= k_cutoff) = 1;
    out = real(ifft(fft(field_in) .* mask));
end

function dfield_dx = spectral_derivative_x_local(field_in, x_vec)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    dfield_dx = real(ifft(1i * kx .* fft(field_in)));
end

function field_out = frequency_filtering_1d_local(field_in, x_vec, kp, n)
    x_vec = x_vec(:);
    field_in = field_in(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    sigma_k = 0.5 * kp;
    k_target = n * kp;
    mask = exp(-((abs(kx) - k_target).^2) / (2 * sigma_k^2));
    field_out = real(ifft(fft(field_in) .* mask));
end

function spec = extract_linear_components_local(eta11, x_vec, CFG)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    nx = numel(x_vec);
    dx = x_vec(2) - x_vec(1);
    dk = 2 * pi / (nx * dx);
    kx_grid = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dk;
    eta_hat = fft(eta11) / nx;

    if mod(nx, 2) == 0
        positive_idx = 2:(nx / 2);
    else
        positive_idx = 2:((nx + 1) / 2);
    end

    ref_amp = max(abs(eta_hat));
    positive_idx = positive_idx(abs(eta_hat(positive_idx)) > CFG.linear_fft_rel_tol * ref_amp);
    if ~CFG.keep_all_positive_modes
        positive_idx = select_energy_dominant_indices_local(eta_hat, positive_idx, CFG.linear_energy_keep);
    end

    coeff = eta_hat(positive_idx);
    spec = struct();
    spec.indices = positive_idx(:);
    spec.k = kx_grid(positive_idx);
    spec.coeff = coeff(:);
    spec.amplitude = 2 * abs(coeff(:));
    spec.phase = angle(coeff(:));
    spec.a = 2 * real(coeff(:));
    spec.b = 2 * imag(coeff(:));
    spec.num_components = numel(positive_idx);
end

function eta1 = reconstruct_first_order_eta_local(spec, x_vec, t_eval)
    x_vec = x_vec(:).';
    eta1 = zeros(size(x_vec));

    for i = 1:spec.num_components
        eta1 = eta1 + spec.amplitude(i) * cos(spec.k(i) * x_vec - 0 * t_eval + spec.phase(i));
    end

    eta1 = eta1(:);
end

function u1 = reconstruct_first_order_u_local(spec, x_vec, z_phys, depth_value, gravity, t_eval)
    x_vec = x_vec(:).';
    z_bed = z_phys + depth_value;
    u1 = zeros(size(x_vec));

    omega = sqrt(gravity * spec.k .* tanh(spec.k * depth_value));
    for i = 1:spec.num_components
        amp = gravity * spec.amplitude(i) * spec.k(i) / omega(i);
        shape = cosh(spec.k(i) * z_bed) / cosh(spec.k(i) * depth_value);
        u1 = u1 + amp * shape * cos(spec.k(i) * x_vec - 0 * t_eval + spec.phase(i));
    end

    u1 = u1(:);
end

function out = compute_strict_difference_frequency_u_local(spec, x_vec, z_phys, depth_value, gravity, t_eval)
    x_vec = x_vec(:).';
    z_bed = z_phys + depth_value;

    k = spec.k(:);
    A = spec.amplitude(:);
    eps = spec.phase(:);
    omega = sqrt(gravity * k .* tanh(k * depth_value));

    u = zeros(size(x_vec));
    pair_meta = zeros(max(spec.num_components * (spec.num_components - 1) / 2, 1), 11);
    pair_count = 0;

    for i = 1:spec.num_components
        for j = (i + 1):spec.num_components
            ki = k(i);
            kj = k(j);
            wi = omega(i);
            wj = omega(j);

            k_minus = ki - kj;
            w_minus = wi - wj;
            D_minus = w_minus^2 - gravity * k_minus * tanh(k_minus * depth_value);

            if abs(D_minus) < 1e-12
                continue;
            end

            A_minus = (wi * wj * w_minus / D_minus) * (1 + 1 / (tanh(ki * depth_value) * tanh(kj * depth_value))) + ...
                (1 / (2 * D_minus)) * (wi^3 / sinh(ki * depth_value)^2 - wj^3 / sinh(kj * depth_value)^2);
            B_minus = (wi^2 + wj^2) / (2 * gravity) + ...
                (wi * wj / (2 * gravity)) * (1 + 1 / (tanh(ki * depth_value) * tanh(kj * depth_value))) * ...
                ((w_minus^2 + gravity * k_minus * tanh(k_minus * depth_value)) / D_minus) + ...
                (w_minus / (2 * gravity * D_minus)) * ...
                (wi^3 / sinh(ki * depth_value)^2 - wj^3 / sinh(kj * depth_value)^2);

            phase = k_minus * x_vec - w_minus * t_eval + (eps(i) - eps(j));
            z_shape = cosh(k_minus * z_bed) / cosh(k_minus * depth_value);
            u = u + A(i) * A(j) * A_minus * k_minus * z_shape .* cos(phase);

            pair_count = pair_count + 1;
            pair_meta(pair_count, :) = [i, j, ki, kj, k_minus, wi, wj, w_minus, D_minus, A_minus, B_minus];
        end
    end

    out = struct();
    out.u = u(:);
    out.num_pairs = pair_count;
    out.k = k;
    out.omega = omega;
    out.amplitude = A;
    out.phase = eps;
    out.z_phys = z_phys;
    out.z_bed = z_bed;
    out.depth = depth_value;
    out.pairs = pair_meta(1:pair_count, :);
    out.pair_columns = {'i', 'j', 'ki', 'kj', 'ki_minus_kj', 'omega_i', 'omega_j', 'omega_minus', 'D_minus', 'A_minus', 'B_minus'};
    out.note = 'Strict Appendix A difference-frequency horizontal velocity only.';
end

function metrics = compute_metrics_local(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);

    if numel(reference) ~= numel(candidate)
        error('Metric input size mismatch.');
    end

    cc = corrcoef(reference, candidate);
    if numel(cc) == 1
        corr_value = 1.0;
    else
        corr_value = cc(1, 2);
    end

    metrics = struct();
    metrics.corr = corr_value;
    metrics.rmse = sqrt(mean((candidate - reference).^2));
    metrics.peak_ratio = max(abs(candidate)) / max(abs(reference));
end

function indices = select_energy_dominant_indices_local(fft_signal, candidate_idx, energy_keep)
    if isempty(candidate_idx)
        indices = candidate_idx;
        return;
    end

    energy = abs(fft_signal(candidate_idx)).^2;
    [energy_sorted, order] = sort(energy, 'descend');
    cumulative = cumsum(energy_sorted) / sum(energy_sorted);
    n_keep = find(cumulative >= energy_keep, 1, 'first');
    keep_idx = sort(candidate_idx(order(1:n_keep)));
    indices = keep_idx(:).';
end

function kd = extract_kd_from_case_pattern_local(folder_pattern)
    token = regexp(folder_pattern, 'kd(?<kd>\d+(?:\.\d+)?)', 'names', 'once');
    if isempty(token)
        error('Unable to parse kd from folder pattern: %s', folder_pattern);
    end
    kd = str2double(token.kd);
end

function label = build_location_label_local(CFG, z_mean)
    switch lower(CFG.sigma_mode)
        case 'surface'
            label = sprintf('surface (mean z = %.4f m)', z_mean);
        case 'z_target'
            label = sprintf('mean z = %.4f m', z_mean);
        otherwise
            label = sprintf('mean z = %.4f m', z_mean);
    end
end

function [x_plot_shifted, eta11_shifted, field1_shifted, field2_shifted, field3_shifted, field4_shifted, field5_shifted, x_center] = ...
        recenter_fields_for_plot_local(x_plot, eta11, eta_phases, field1, field2, field3, field4, field5, window_source)
    x_plot = x_plot(:);
    eta11 = eta11(:);
    field1 = field1(:);
    field2 = field2(:);
    field3 = field3(:);
    field4 = field4(:);
    field5 = field5(:);

    switch lower(window_source)
        case 'phase0_eta'
            window_field = eta_phases(1, :).';
        case 'eta11'
            window_field = eta11;
        otherwise
            error('Unsupported CFG.window_source: %s. Use ''phase0_eta'' or ''eta11''.', window_source);
    end

    if max(abs(window_field)) < 1e-12
        window_field = field1;
    end

    shift_idx = compute_window_shift_local(window_field);
    x_plot_shifted = x_plot;
    eta11_shifted = circshift(eta11, shift_idx);
    field1_shifted = circshift(field1, shift_idx);
    field2_shifted = circshift(field2, shift_idx);
    field3_shifted = circshift(field3, shift_idx);
    field4_shifted = circshift(field4, shift_idx);
    field5_shifted = circshift(field5, shift_idx);

    recentered_window_field = circshift(window_field(:), shift_idx);
    recentered_envelope = abs(hilbert(recentered_window_field));
    [~, recentered_peak_idx] = max(recentered_envelope);
    x_center = x_plot_shifted(recentered_peak_idx);
end

function field_shifted = shift_with_window_local(field_in, eta_phases, eta11, window_source)
    field_in = field_in(:);

    switch lower(window_source)
        case 'phase0_eta'
            window_field = eta_phases(1, :).';
        case 'eta11'
            window_field = eta11(:);
        otherwise
            error('Unsupported CFG.window_source: %s.', window_source);
    end

    shift_idx = compute_window_shift_local(window_field);
    field_shifted = circshift(field_in, shift_idx);
end

function shift_idx = compute_window_shift_local(window_field)
    window_envelope = abs(hilbert(window_field(:)));
    [~, peak_idx] = max(window_envelope);
    n = numel(window_field);
    center_idx = floor((n + 1) / 2);
    shift_idx = center_idx - peak_idx;
end

function [k_pos, amp_a, amp_b, amp_c] = compute_single_sided_spectrum_local(field_a, field_b, field_c, x_vec)
    field_a = field_a(:);
    field_b = field_b(:);
    field_c = field_c(:);
    x_vec = x_vec(:);

    nx = numel(x_vec);
    dx = x_vec(2) - x_vec(1);
    dk = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dk;

    if mod(nx, 2) == 0
        pos_idx = 1:(nx / 2 + 1);
    else
        pos_idx = 1:((nx + 1) / 2);
    end

    fa = fft(field_a) / nx;
    fb = fft(field_b) / nx;
    fc = fft(field_c) / nx;

    amp_a = 2 * abs(fa(pos_idx));
    amp_b = 2 * abs(fb(pos_idx));
    amp_c = 2 * abs(fc(pos_idx));
    amp_a(1) = abs(fa(1));
    amp_b(1) = abs(fb(1));
    amp_c(1) = abs(fc(1));

    if mod(nx, 2) == 0
        amp_a(end) = abs(fa(pos_idx(end)));
        amp_b(end) = abs(fb(pos_idx(end)));
        amp_c(end) = abs(fc(pos_idx(end)));
    end

    k_pos = abs(kx(pos_idx));
end
