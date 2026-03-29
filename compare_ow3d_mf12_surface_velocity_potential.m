% compare_ow3d_mf12_surface_velocity_potential.m
% Compare OW3D and MF12 free-surface velocity potential using the existing
% four-phase harmonic-separation workflow.

clc;
clear;
close all;

addpath(fullfile(pwd, 'irregularWavesMF12', 'Source'));
addpath(fullfile(pwd, 'test functions for VWA Opensource'));

CFG = struct();
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check3');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.kinematics_file_id = 1;
CFG.time_index = [];
CFG.default_time_index_from_end = 180;
CFG.lambda = 225;
CFG.gravity = 9.81;
CFG.kp_depth = 0.0279;
CFG.apply_x_filter = false;
CFG.apply_eta11_bandpass = false;
CFG.linear_fft_rel_tol = 1e-12;
CFG.keep_all_positive_modes = false;
CFG.linear_energy_keep = 0.999;
CFG.linear_min_components = 0;
CFG.max_linear_components = 380;
CFG.plot_window_lambda = 5.0;
CFG.filter_band_subharmonic = [0.0, 1.5];
CFG.filter_band_order1 = [0.0, 3.0];
CFG.filter_band_order2 = [0.8, 3.5];
CFG.filter_band_order3 = [1.5, 5.0];
CFG.filter_transition_low = 0.25;
CFG.filter_transition_high = 0.35;
CFG.subharmonic_cutoff_factor = 0.75;
CFG.subharmonic_transition_factor = 1.00;
CFG.include_vwa = true;
CFG.vwa_small_kd_min = 0.3;
CFG.vwa_analytic_side = 'neg';
CFG.vwa_phi_take = 'imag';
CFG.mf12_disable_third_order_correction = true;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');
CFG.save_mat = false;

four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

phase_data = cell(1, numel(CFG.phi_shifts_deg));
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
surface_sigma_idx = numel(sigma_vec);
t_selected = phase_data{1}.t(selected_time_index);
kp = 2 * pi / CFG.lambda;
case_kd = extract_kd_from_case_pattern_local(CFG.folder_pattern);
depth_value = case_kd / CFG.kp_depth;

eta_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));
phi_surface_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec));
for idx = 1:numel(CFG.phi_shifts_deg)
    eta_phases(idx, :) = squeeze(phase_data{idx}.eta(selected_time_index, :, 1));
    phi_surface_phases(idx, :) = squeeze(phase_data{idx}.phi(selected_time_index, surface_sigma_idx, :, 1));
end

eta_harmonics_raw = reconstruct_harmonics_1d_hilbert_local(eta_phases, four_phase_coef);
phi_harmonics_raw = reconstruct_harmonics_1d_hilbert_local(phi_surface_phases, four_phase_coef);

eta_harmonics = eta_harmonics_raw;
phi_harmonics = phi_harmonics_raw;

if CFG.apply_x_filter
    eta_harmonics = filter_harmonics_x_only_local(eta_harmonics, x_vec, kp);
    phi_harmonics = filter_harmonics_x_only_local(phi_harmonics, x_vec, kp);
end

eta11 = eta_harmonics(1, :).';
if CFG.apply_eta11_bandpass
    eta11 = frequency_filtering_1d_local(eta11, x_vec, kp, 1);
end

linear_spec = extract_linear_components_local(eta11, x_vec, CFG);
print_linear_component_summary_local(linear_spec, kp);
mf12_surface_raw = compute_mf12_surface_phi_orders_local(linear_spec, x_vec, depth_value, CFG.gravity, CFG);
mf12_surface = mf12_surface_raw;
if CFG.apply_x_filter
    mf12_surface = filter_surface_orders_x_only_local(mf12_surface_raw, x_vec, kp);
end

vwa_surface_raw = struct();
vwa_surface = struct();
vwa_components_raw = struct();
vwa_components = struct();
metrics_vwa = struct();
metrics_mf12_vwa = struct();

if CFG.include_vwa
    vwa_surface_raw = compute_vwa_surface_phi_orders_local(eta11, x_vec, depth_value, CFG.gravity, CFG);
    vwa_surface = vwa_surface_raw;
    vwa_components_raw = build_surface_comparison_components_local(vwa_surface_raw, x_vec, kp, CFG, 'vwa');
    vwa_components = build_surface_comparison_components_local(vwa_surface, x_vec, kp, CFG, 'vwa');
end

eta_metrics = struct();
eta_metrics.order1 = compute_metrics_local(eta11, mf12_surface_raw.eta1);
eta_metrics.order1_filtered = compute_metrics_local(eta11, mf12_surface.eta1);

ow3d_surface_raw = struct();
ow3d_surface_raw.order1 = phi_harmonics_raw(1, :).';
ow3d_surface_raw.order2 = phi_harmonics_raw(2, :).';
ow3d_surface_raw.order3 = phi_harmonics_raw(3, :).';
ow3d_surface_raw.order4 = phi_harmonics_raw(4, :).';

ow3d_surface = struct();
ow3d_surface.order1 = phi_harmonics(1, :).';
ow3d_surface.order2 = phi_harmonics(2, :).';
ow3d_surface.order3 = phi_harmonics(3, :).';
ow3d_surface.order4 = phi_harmonics(4, :).';

ow3d_components_raw = build_surface_comparison_components_local(ow3d_surface_raw, x_vec, kp, CFG, 'ow3d');
mf12_components_raw = build_surface_comparison_components_local(mf12_surface_raw, x_vec, kp, CFG, 'mf12');

ow3d_components = build_surface_comparison_components_local(ow3d_surface, x_vec, kp, CFG, 'ow3d');
mf12_components = build_surface_comparison_components_local(mf12_surface, x_vec, kp, CFG, 'mf12');

metrics = struct();
metrics.first = compute_metrics_local(ow3d_components.first, mf12_components.first);
metrics.second_super = compute_metrics_local(ow3d_components.second_super, mf12_components.second_super);
metrics.third = compute_metrics_local(ow3d_components.third, mf12_components.third);
metrics.second_sub = compute_metrics_local(ow3d_components.second_sub, mf12_components.second_sub);
if CFG.include_vwa
    metrics_vwa.second_super = compute_metrics_local(ow3d_components.second_super, vwa_components.second_super);
    metrics_vwa.third = compute_metrics_local(ow3d_components.third, vwa_components.third);
    metrics_mf12_vwa.second_super = compute_metrics_local(mf12_components.second_super, vwa_components.second_super);
    metrics_mf12_vwa.third = compute_metrics_local(mf12_components.third, vwa_components.third);
end

surface_split = compute_surface_frequency_split_local( ...
    ow3d_components_raw, mf12_components_raw, x_vec, kp, ...
    CFG.subharmonic_cutoff_factor * kp, CFG.subharmonic_transition_factor * kp);

x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / CFG.lambda;
[x_plot_shifted, eta11_shifted, phi1_ow3d_shifted, phi1_mf12_shifted, phi2sup_ow3d_shifted, phi2sup_mf12_shifted, x_center_window] = ...
    recenter_fields_for_plot_local(x_plot, eta11, eta_phases, ow3d_components.first, mf12_components.first, ...
    ow3d_components.second_super, mf12_components.second_super);
eta1_mf12_shifted = shift_with_window_local(mf12_surface.eta1, eta_phases, eta11);
phi3_ow3d_shifted = shift_with_window_local(ow3d_components.third, eta_phases, eta11);
phi3_mf12_shifted = shift_with_window_local(mf12_components.third, eta_phases, eta11);
phi2sub_ow3d_shifted = shift_with_window_local(ow3d_components.second_sub, eta_phases, eta11);
phi2sub_mf12_shifted = shift_with_window_local(mf12_components.second_sub, eta_phases, eta11);
phi2sup_vwa_shifted = [];
phi3_vwa_shifted = [];
if CFG.include_vwa
    phi2sup_vwa_shifted = shift_with_window_local(vwa_components.second_super, eta_phases, eta11);
    phi3_vwa_shifted = shift_with_window_local(vwa_components.third, eta_phases, eta11);
end

x_limits = [x_center_window - 0.5 * CFG.plot_window_lambda, x_center_window + 0.5 * CFG.plot_window_lambda];
x_limits(1) = max(x_limits(1), x_plot_shifted(1));
x_limits(2) = min(x_limits(2), x_plot_shifted(end));

y_limits = compute_pairwise_ylimits_local( ...
    {phi1_ow3d_shifted, phi2sup_ow3d_shifted, phi3_ow3d_shifted, phi2sub_ow3d_shifted}, ...
    {phi1_mf12_shifted, phi2sup_mf12_shifted, phi3_mf12_shifted, phi2sub_mf12_shifted});

if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

fig = figure('Color', 'w', 'Position', [120 60 1500 1080]);
tile = tiledlayout(fig, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile, sprintf('Surface velocity potential: OW3D vs MF12 (t = %.4f s, kd = %.2f)', t_selected, case_kd), ...
    'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

series_labels = {'First harmonic', 'Second superharmonic', 'Third harmonic', 'Second subharmonic'};
ow3d_series = {phi1_ow3d_shifted, phi2sup_ow3d_shifted, phi3_ow3d_shifted, phi2sub_ow3d_shifted};
mf12_series = {phi1_mf12_shifted, phi2sup_mf12_shifted, phi3_mf12_shifted, phi2sub_mf12_shifted};
metric_series = {metrics.first, metrics.second_super, metrics.third, metrics.second_sub};
panel_prefix = {'(a)', '(b)', '(c)', '(d)'};

for panel_idx = 1:4
    ax = nexttile(tile);
    draw_comparison_panel_local(ax, x_plot_shifted, ow3d_series{panel_idx}, mf12_series{panel_idx}, ...
        sprintf('%s %s', panel_prefix{panel_idx}, series_labels{panel_idx}), ...
        y_limits(panel_idx, :), '$\widetilde{\Phi}$ (m$^2$/s)', x_limits);
    text(ax, 0.02, 0.92, sprintf('corr = %.3f, RMSE = %.3e, peak ratio = %.3f', ...
        metric_series{panel_idx}.corr, metric_series{panel_idx}.rmse, metric_series{panel_idx}.peak_ratio), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);
end

annotation(fig, 'textbox', [0.13 0.01 0.82 0.04], ...
    'String', sprintf(['OW3D reference: surface sigma layer from four-phase kinematics reconstruction. ' ...
    'MF12 input: linear spectrum extracted from OW3D \\eta^{(1)} and evaluated with mf12\\_spectral\\_surface. ' ...
    'Displayed second-order content is split into superharmonic and subharmonic bands. ' ...
    'Retained %d of %d linear components for MF12 reconstruction.'], ...
    mf12_surface_raw.spec_used.num_components, linear_spec.num_components), ...
    'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'FontName', 'Times New Roman', 'FontSize', 11);

png_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_surface_velocity_potential.png');
exportgraphics(fig, png_path, 'Resolution', 300);

png_super_path = '';
png_super_spec_path = '';
if CFG.include_vwa
    y_limits_super = compute_threeway_ylimits_local( ...
        {phi2sup_ow3d_shifted, phi3_ow3d_shifted}, ...
        {phi2sup_mf12_shifted, phi3_mf12_shifted}, ...
        {phi2sup_vwa_shifted, phi3_vwa_shifted});

    fig_super = figure('Color', 'w', 'Position', [180 80 1500 760]);
    tile_super = tiledlayout(fig_super, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile_super, sprintf('Superharmonic surface velocity potential: OW3D vs MF12 vs VWA (t = %.4f s, kd = %.2f)', ...
        t_selected, case_kd), 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    ax = nexttile(tile_super);
    draw_threeway_comparison_panel_local(ax, x_plot_shifted, phi2sup_ow3d_shifted, phi2sup_mf12_shifted, phi2sup_vwa_shifted, ...
        '(a) Second superharmonic', y_limits_super(1, :), '$\widetilde{\Phi}$ (m$^2$/s)', x_limits, ...
        'OW3D', 'MF12', 'VWA');
    text(ax, 0.02, 0.92, sprintf('OW3D-MF12 corr = %.3f | OW3D-VWA corr = %.3f | MF12-VWA corr = %.3f', ...
        metrics.second_super.corr, metrics_vwa.second_super.corr, metrics_mf12_vwa.second_super.corr), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

    ax = nexttile(tile_super);
    draw_threeway_comparison_panel_local(ax, x_plot_shifted, phi3_ow3d_shifted, phi3_mf12_shifted, phi3_vwa_shifted, ...
        '(b) Third harmonic', y_limits_super(2, :), '$\widetilde{\Phi}$ (m$^2$/s)', x_limits, ...
        'OW3D', 'MF12', 'VWA');
    text(ax, 0.02, 0.92, sprintf('OW3D-MF12 corr = %.3f | OW3D-VWA corr = %.3f | MF12-VWA corr = %.3f', ...
        metrics.third.corr, metrics_vwa.third.corr, metrics_mf12_vwa.third.corr), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

    annotation(fig_super, 'textbox', [0.13 0.01 0.82 0.04], ...
        'String', sprintf(['VWA input uses the same OW3D \\eta^{(1)} as MF12. ' ...
        'VWA surface potential is evaluated with vwa\\_compute orders 2-3 only, so this three-way figure is restricted to superharmonic content.']), ...
        'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontName', 'Times New Roman', 'FontSize', 11);

    png_super_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_vwa_surface_velocity_potential_superharmonics.png');
    exportgraphics(fig_super, png_super_path, 'Resolution', 300);
end

fig_eta = figure('Color', 'w', 'Position', [160 120 1500 620]);
tile_eta = tiledlayout(fig_eta, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile_eta, sprintf('First-harmonic free-surface elevation: OW3D vs MF12 (t = %.4f s, kd = %.2f)', ...
    t_selected, case_kd), 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

eta_y_limits = compute_pairwise_ylimits_local({eta11_shifted}, {eta1_mf12_shifted});
ax_eta = nexttile(tile_eta);
draw_comparison_panel_local(ax_eta, x_plot_shifted, eta11_shifted, eta1_mf12_shifted, ...
    '(a) First-harmonic free-surface elevation', eta_y_limits(1, :), '\eta^{(1)} (m)', x_limits);
text(ax_eta, 0.02, 0.92, sprintf('corr = %.3f, RMSE = %.3e, peak ratio = %.3f', ...
    eta_metrics.order1_filtered.corr, eta_metrics.order1_filtered.rmse, eta_metrics.order1_filtered.peak_ratio), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

ax_eta_diff = nexttile(tile_eta);
eta_diff = eta11_shifted - eta1_mf12_shifted;
eta_diff_limits = compute_pairwise_ylimits_local({eta_diff}, {zeros(size(eta_diff))});
draw_difference_panel_local(ax_eta_diff, x_plot_shifted, eta_diff, ...
    '(b) Difference: OW3D - MF12', eta_diff_limits(1, :), '\Delta \eta^{(1)} (m)', x_limits);

png_eta_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_eta1.png');
exportgraphics(fig_eta, png_eta_path, 'Resolution', 300);

[k_plot, amp_phi1_ow3d, amp_phi1_mf12] = compute_spectrum_pair_local( ...
    ow3d_components.first, mf12_components.first, x_vec);
[~, amp_phi2sup_ow3d, amp_phi2sup_mf12] = compute_spectrum_pair_local( ...
    ow3d_components.second_super, mf12_components.second_super, x_vec);
[~, amp_phi3_ow3d, amp_phi3_mf12] = compute_spectrum_pair_local( ...
    ow3d_components.third, mf12_components.third, x_vec);
[~, amp_phi2sub_ow3d, amp_phi2sub_mf12] = compute_spectrum_pair_local( ...
    ow3d_components.second_sub, mf12_components.second_sub, x_vec);
amp_phi2sup_vwa = [];
amp_phi3_vwa = [];
if CFG.include_vwa
    [~, ~, ~, amp_phi2sup_vwa, amp_phi3_vwa, ~] = compute_three_pair_spectra_local( ...
        ow3d_components.second_super, mf12_components.second_super, vwa_components.second_super, ...
        ow3d_components.third, mf12_components.third, vwa_components.third, x_vec);
end

fig_spec = figure('Color', 'w', 'Position', [150 80 1500 1080]);
tile = tiledlayout(fig_spec, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile, sprintf('Surface velocity potential spectra: OW3D vs MF12 (t = %.4f s)', t_selected), ...
    'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

plot_spectrum_pair_local(nexttile(tile), k_plot / kp, amp_phi1_ow3d, amp_phi1_mf12, '(a) First-harmonic spectrum');
plot_spectrum_pair_local(nexttile(tile), k_plot / kp, amp_phi2sup_ow3d, amp_phi2sup_mf12, '(b) Second-superharmonic spectrum');
plot_spectrum_pair_local(nexttile(tile), k_plot / kp, amp_phi3_ow3d, amp_phi3_mf12, '(c) Third-harmonic spectrum');
plot_spectrum_pair_local(nexttile(tile), k_plot / kp, amp_phi2sub_ow3d, amp_phi2sub_mf12, '(d) Second-subharmonic spectrum');

png_spec_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_surface_velocity_potential_spectra.png');
exportgraphics(fig_spec, png_spec_path, 'Resolution', 300);

if CFG.include_vwa
    fig_super_spec = figure('Color', 'w', 'Position', [220 100 1500 760]);
    tile_super_spec = tiledlayout(fig_super_spec, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile_super_spec, sprintf('Superharmonic surface-potential spectra: OW3D vs MF12 vs VWA (t = %.4f s)', ...
        t_selected), 'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    plot_spectrum_threeway_local(nexttile(tile_super_spec), k_plot / kp, amp_phi2sup_ow3d, amp_phi2sup_mf12, amp_phi2sup_vwa, ...
        '(a) Second-superharmonic spectrum');
    plot_spectrum_threeway_local(nexttile(tile_super_spec), k_plot / kp, amp_phi3_ow3d, amp_phi3_mf12, amp_phi3_vwa, ...
        '(b) Third-harmonic spectrum');

    png_super_spec_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_vwa_surface_velocity_potential_superharmonic_spectra.png');
    exportgraphics(fig_super_spec, png_super_spec_path, 'Resolution', 300);
end

results = struct();
results.cfg = CFG;
results.selected_time_index = selected_time_index;
results.t_selected = t_selected;
results.kp = kp;
results.kd = case_kd;
results.depth_value = depth_value;
results.x = x_vec;
results.x_plot = x_plot;
results.surface_sigma_idx = surface_sigma_idx;
results.surface_sigma_value = sigma_vec(surface_sigma_idx);
results.eta11 = eta11;
results.linear_spec = linear_spec;
results.eta_metrics = eta_metrics;
results.ow3d_surface_raw = ow3d_surface_raw;
results.ow3d_surface = ow3d_surface;
results.mf12_surface_raw = mf12_surface_raw;
results.mf12_surface = mf12_surface;
results.vwa_surface_raw = vwa_surface_raw;
results.vwa_surface = vwa_surface;
results.ow3d_components = ow3d_components;
results.ow3d_components_raw = ow3d_components_raw;
results.mf12_components = mf12_components;
results.mf12_components_raw = mf12_components_raw;
results.vwa_components = vwa_components;
results.vwa_components_raw = vwa_components_raw;
results.metrics = metrics;
results.metrics_vwa = metrics_vwa;
results.metrics_mf12_vwa = metrics_mf12_vwa;
results.frequency_split = surface_split;
results.spectra = struct( ...
    'k', k_plot, ...
    'phi1_ow3d', amp_phi1_ow3d, 'phi1_mf12', amp_phi1_mf12, ...
    'phi2super_ow3d', amp_phi2sup_ow3d, 'phi2super_mf12', amp_phi2sup_mf12, ...
    'phi2super_vwa', amp_phi2sup_vwa, ...
    'phi3_ow3d', amp_phi3_ow3d, 'phi3_mf12', amp_phi3_mf12, ...
    'phi3_vwa', amp_phi3_vwa, ...
    'phi2sub_ow3d', amp_phi2sub_ow3d, 'phi2sub_mf12', amp_phi2sub_mf12);
results.eta1_compare = struct( ...
    'ow3d', eta11(:), ...
    'mf12', mf12_surface.eta1(:), ...
    'ow3d_shifted', eta11_shifted(:), ...
    'mf12_shifted', eta1_mf12_shifted(:));

mat_path = '';
if CFG.save_mat
    mat_path = fullfile(CFG.output_dir, 'compare_ow3d_mf12_surface_velocity_potential.mat');
    save(mat_path, 'results');
end

fprintf('\n=== Surface velocity potential: OW3D vs MF12 ===\n');
fprintf('Eta(1) raw MF12 vs input eta11: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    eta_metrics.order1.corr, eta_metrics.order1.rmse, eta_metrics.order1.peak_ratio);
fprintf('Eta(1) filtered MF12 vs input eta11: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    eta_metrics.order1_filtered.corr, eta_metrics.order1_filtered.rmse, eta_metrics.order1_filtered.peak_ratio);
fprintf('First harmonic: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    metrics.first.corr, metrics.first.rmse, metrics.first.peak_ratio);
fprintf('Second superharmonic: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    metrics.second_super.corr, metrics.second_super.rmse, metrics.second_super.peak_ratio);
fprintf('Third harmonic: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    metrics.third.corr, metrics.third.rmse, metrics.third.peak_ratio);
fprintf('Second subharmonic: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
    metrics.second_sub.corr, metrics.second_sub.rmse, metrics.second_sub.peak_ratio);
if CFG.include_vwa
    fprintf('\n=== Superharmonic comparison including VWA ===\n');
    fprintf('Second superharmonic OW3D vs VWA: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
        metrics_vwa.second_super.corr, metrics_vwa.second_super.rmse, metrics_vwa.second_super.peak_ratio);
    fprintf('Second superharmonic MF12 vs VWA: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
        metrics_mf12_vwa.second_super.corr, metrics_mf12_vwa.second_super.rmse, metrics_mf12_vwa.second_super.peak_ratio);
    fprintf('Third harmonic OW3D vs VWA: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
        metrics_vwa.third.corr, metrics_vwa.third.rmse, metrics_vwa.third.peak_ratio);
    fprintf('Third harmonic MF12 vs VWA: corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
        metrics_mf12_vwa.third.corr, metrics_mf12_vwa.third.rmse, metrics_mf12_vwa.third.peak_ratio);
end
fprintf('\n=== Frequency-separated comparison using raw surface potential ===\n');
fprintf('Lowpass cutoff: %.3f kp, transition: %.3f kp\n', ...
    CFG.subharmonic_cutoff_factor, CFG.subharmonic_transition_factor);
print_frequency_split_metrics_local(surface_split);
fprintf('Saved figure: %s\n', png_path);
if ~isempty(png_super_path)
    fprintf('Saved superharmonic three-way figure: %s\n', png_super_path);
end
fprintf('Saved eta(1) figure: %s\n', png_eta_path);
fprintf('Saved spectra figure: %s\n', png_spec_path);
if ~isempty(png_super_spec_path)
    fprintf('Saved superharmonic spectra figure: %s\n', png_super_spec_path);
end
if ~isempty(mat_path)
    fprintf('Saved MAT: %s\n', mat_path);
end

function kin_path = resolve_kinematics_path_local(case_folder, file_id)
    if file_id < 10
        file_name = sprintf('Kinematics0%d.bin', file_id);
    else
        file_name = sprintf('Kinematics%d.bin', file_id);
    end
    kin_path = fullfile(case_folder, file_name);
end

function time_index = resolve_time_index_local(n_times_valid, cfg_time_index, default_time_index_from_end)
    if isempty(cfg_time_index)
        time_index = n_times_valid - default_time_index_from_end + 1;
    elseif isscalar(cfg_time_index) && cfg_time_index < 0
        time_index = n_times_valid + cfg_time_index + 1;
    else
        time_index = cfg_time_index;
    end
    time_index = round(time_index(1));
    time_index = min(max(1, time_index), n_times_valid);
end

function data = read_ow3d_kinematics_full_local(file_path)
    [it, eta, phi, x, y, h, sigma, t] = read_kinematics_file_core_local(file_path);

    data = struct();
    data.it = it;
    data.n_times_valid = it;
    data.eta = eta;
    data.phi = phi;
    data.x = x;
    data.y = y;
    data.h = h;
    data.sigma = sigma;
    data.t = t;
end

function [it, eta, phi, x, y, h, sigma, t] = read_kinematics_file_core_local(file_path)
    fid = fopen(file_path, 'r', 'ieee-le');
    if fid < 0
        error('Could not open kinematics file: %s', file_path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fread(fid, 1, 'int');
    xbeg = fread(fid, 1, 'int');
    xend = fread(fid, 1, 'int');
    xstride = fread(fid, 1, 'int');
    ybeg = fread(fid, 1, 'int');
    yend = fread(fid, 1, 'int');
    ystride = fread(fid, 1, 'int');
    tbeg = fread(fid, 1, 'int');
    tend = fread(fid, 1, 'int');
    tstride = fread(fid, 1, 'int');
    dt = fread(fid, 1, 'double');
    nz = fread(fid, 1, 'int');
    sigma = zeros(nz, 1);
    fread(fid, 2, 'int');

    nx = floor((xend - xbeg) / xstride) + 1;
    ny = floor((yend - ybeg) / ystride) + 1;
    nt = floor((tend - tbeg) / tstride) + 1;

    header = fread(fid, 5 * nx * ny, 'double');
    fread(fid, 2, 'int');

    x = zeros(nx, ny);
    y = zeros(nx, ny);
    h = zeros(nx, ny);
    x(:) = header(1:5:end);
    y(:) = header(2:5:end);
    h(:) = header(3:5:end);

    for i = 1:nz
        sigma(i) = fread(fid, 1, 'double');
    end
    fread(fid, 2, 'int');

    eta = zeros(nt, nx, ny);
    phi = zeros(nt, nz, nx, ny);
    t = (0:nt - 1) * dt * tstride;
    it = 0;

    for it_read = 1:nt - 1
        tmp_eta = fread(fid, nx * ny, 'double');
        if numel(tmp_eta) < nx * ny
            it = it_read - 1;
            break;
        end
        eta(it_read, :) = tmp_eta;
        fread(fid, 2, 'int');

        tmp_skip = fread(fid, nx * ny, 'double');
        if numel(tmp_skip) < nx * ny
            it = it_read - 1;
            break;
        end
        fread(fid, 2, 'int');

        tmp_skip = fread(fid, nx * ny, 'double');
        if numel(tmp_skip) < nx * ny
            it = it_read - 1;
            break;
        end
        fread(fid, 2, 'int');

        tmp_phi = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_phi) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        phi(it_read, :) = tmp_phi;
        fread(fid, 2, 'int');

        for block_idx = 1:6
            tmp_skip = fread(fid, nx * ny * nz, 'double');
            if numel(tmp_skip) < nx * ny * nz
                it = it_read - 1;
                break;
            end
            fread(fid, 2, 'int');
        end
        if numel(tmp_skip) < nx * ny * nz
            break;
        end
        it = it_read;
    end

    if it <= 0
        error('No complete stored kinematics time step could be read from %s', file_path);
    end

    eta = eta(1:it, :, :);
    phi = phi(1:it, :, :, :);
    t = t(1:it);
end

function harmonics = reconstruct_harmonics_1d_hilbert_local(fields_by_phase, coef)
    analytic_part = hilbert(fields_by_phase.').';
    all_fields = cat(1, real(fields_by_phase), -imag(analytic_part));
    harmonics = zeros(4, size(fields_by_phase, 2));
    for n = 1:4
        harmonics(n, :) = coef(n, :) * all_fields;
    end
end

function harmonics_out = filter_harmonics_x_only_local(harmonics_in, x_vec, kp)
    harmonics_out = harmonics_in;
    for n = 1:size(harmonics_in, 1)
        harmonics_out(n, :) = frequency_filtering_1d_local(squeeze(harmonics_in(n, :)), x_vec, kp, n);
    end
end

function field_out = frequency_filtering_1d_local(field_in, x_vec, kp, n)
    x_vec = x_vec(:);
    field_in = field_in(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    [kmin_factor, kmax_factor] = harmonic_filter_band_local(n);
    transition_low = 0.25 * kp;
    transition_high = 0.35 * kp;
    mask = smooth_bandpass_mask_local(abs(kx), kmin_factor * kp, kmax_factor * kp, transition_low, transition_high);
    field_out = real(ifft(fft(field_in) .* mask));
end

function [kmin_factor, kmax_factor] = harmonic_filter_band_local(order)
    switch order
        case 0
            kmin_factor = 0.0;
            kmax_factor = 2.0;
        case 1
            kmin_factor = 0.0;
            kmax_factor = 3.0;
        case 2
            kmin_factor = 0.8;
            kmax_factor = 3.5;
        case 3
            kmin_factor = 1.5;
            kmax_factor = 5.0;
        otherwise
            kmin_factor = max(0.0, order - 0.75);
            kmax_factor = order + 2.0;
    end
end

function mask = smooth_bandpass_mask_local(k_abs, kmin, kmax, transition_low, transition_high)
    k_abs = k_abs(:);
    mask = zeros(size(k_abs));

    if kmax <= 0
        return;
    end

    if kmin <= 0
        mask(k_abs <= kmax) = 1.0;
    else
        mask(k_abs >= kmin & k_abs <= kmax) = 1.0;
        if transition_low > 0
            idx = k_abs > max(0, kmin - transition_low) & k_abs < kmin;
            xi = (k_abs(idx) - (kmin - transition_low)) / transition_low;
            mask(idx) = 0.5 - 0.5 * cos(pi * xi);
        end
    end

    if transition_high > 0
        idx = k_abs > kmax & k_abs < (kmax + transition_high);
        xi = (k_abs(idx) - kmax) / transition_high;
        mask(idx) = 0.5 + 0.5 * cos(pi * xi);
    end
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
    candidate_idx = positive_idx(abs(eta_hat(positive_idx)) > CFG.linear_fft_rel_tol * ref_amp);
    positive_idx = enforce_min_component_count_local(eta_hat, positive_idx, candidate_idx, CFG.linear_min_components);
    if ~CFG.keep_all_positive_modes
        positive_idx = select_energy_dominant_indices_local(eta_hat, positive_idx, CFG.linear_energy_keep);
        positive_idx = enforce_min_component_count_local(eta_hat, candidate_idx, positive_idx, CFG.linear_min_components);
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

function print_linear_component_summary_local(spec, kp)
    fprintf('\n=== Linear component retention summary ===\n');
    fprintf('Retained positive-k components: %d\n', spec.num_components);

    if spec.num_components == 0
        fprintf('No linear components were retained.\n');
        return;
    end

    k_norm = spec.k(:) / kp;
    amp = spec.amplitude(:);

    fprintf('k/kp range: [%.6f, %.6f]\n', min(k_norm), max(k_norm));
    fprintf('Amplitude range: [%.6e, %.6e]\n', min(amp), max(amp));
    fprintf('Top retained components by amplitude:\n');

    [amp_sorted, order] = sort(amp, 'descend');
    n_show = min(12, numel(order));
    for i = 1:n_show
        idx = order(i);
        fprintf('  %2d: k/kp = %.6f, amp = %.6e, phase = %.6f rad\n', ...
            i, k_norm(idx), amp_sorted(i), spec.phase(idx));
    end
end

function kept_idx = enforce_min_component_count_local(eta_hat, full_positive_idx, current_idx, min_components)
    if nargin < 4 || isempty(min_components)
        min_components = 0;
    end

    min_components = max(0, round(min_components));
    if min_components == 0
        kept_idx = current_idx;
        return;
    end

    if numel(current_idx) >= min_components
        kept_idx = current_idx;
        return;
    end

    if isempty(full_positive_idx)
        kept_idx = current_idx;
        return;
    end

    spectral_energy = abs(eta_hat(full_positive_idx)).^2;
    [~, order] = sort(spectral_energy, 'descend');
    n_keep = min(min_components, numel(full_positive_idx));
    kept_idx = sort(full_positive_idx(order(1:n_keep)));
end

function out = compute_mf12_surface_phi_orders_local(spec, x_vec, depth_value, gravity, CFG)
    x_vec = x_vec(:);
    t_eval = 0.0;
    nx = numel(x_vec);

    out = struct();
    if isempty(spec.k)
        zeros_field = zeros(nx, 1);
        out.eta1 = zeros_field;
        out.eta2 = zeros_field;
        out.eta3 = zeros_field;
        out.eta_nonlinear = zeros_field;
        out.order1 = zeros_field;
        out.order2 = zeros_field;
        out.order3 = zeros_field;
        out.nonlinear = zeros_field;
        return;
    end

    dx = x_vec(2) - x_vec(1);
    Lx = dx * nx;
    Ly = 1.0;
    Ny = 1;

    spec_used = limit_linear_components_local(spec, CFG.max_linear_components);

    fprintf('MF12 component usage: retained %d of %d linear components\n', ...
        spec_used.num_components, spec.num_components);

    kx1 = spec_used.k(:).';
    ky1 = zeros(size(kx1));
    a1 = spec_used.a(:).';
    b1 = spec_used.b(:).';
    coeffs1 = mf12_spectral_coefficients(1, gravity, depth_value, a1, b1, kx1, ky1, 0, 0);
    [eta1_grid, phi1_grid] = mf12_spectral_surface(coeffs1, Lx, Ly, nx, Ny, t_eval);

    kx2 = spec_used.k(:).';
    ky2 = zeros(size(kx2));
    a2 = spec_used.a(:).';
    b2 = spec_used.b(:).';
    coeffs2 = mf12_spectral_coefficients(2, gravity, depth_value, a2, b2, kx2, ky2, 0, 0);
    [eta12_grid, phi12_grid] = mf12_spectral_surface(coeffs2, Lx, Ly, nx, Ny, t_eval);

    kx3 = spec_used.k(:).';
    ky3 = zeros(size(kx3));
    a3 = spec_used.a(:).';
    b3 = spec_used.b(:).';
    coeffs3 = mf12_spectral_coefficients(3, gravity, depth_value, a3, b3, kx3, ky3, 0, 0, ...
        struct('disable_third_order_correction', CFG.mf12_disable_third_order_correction));
    [eta123_grid, phi123_grid] = mf12_spectral_surface(coeffs3, Lx, Ly, nx, Ny, t_eval);

    eta1 = reshape(eta1_grid.', [], 1);
    eta12 = reshape(eta12_grid.', [], 1);
    eta123 = reshape(eta123_grid.', [], 1);
    phi1 = reshape(phi1_grid.', [], 1);
    phi12 = reshape(phi12_grid.', [], 1);
    phi123 = reshape(phi123_grid.', [], 1);

    out.eta1 = eta1(:);
    out.eta2 = eta12(:) - eta1(:);
    out.eta3 = eta123(:) - eta12(:);
    out.eta_nonlinear = eta123(:);
    out.order1 = phi1(:);
    out.order2 = phi12(:) - phi1(:);
    out.order3 = phi123(:) - phi12(:);
    out.nonlinear = phi123(:);
    out.phi12 = phi12(:);
    out.coeffs1 = coeffs1;
    out.coeffs2 = coeffs2;
    out.coeffs3 = coeffs3;
    out.spec_used = spec_used;
end

function out = compute_vwa_surface_phi_orders_local(eta11, x_vec, depth_value, gravity, CFG)
    vwa_out = vwa_compute(eta11, x_vec, depth_value, gravity, ...
        'nList', [2 3], ...
        'analytic_side', CFG.vwa_analytic_side, ...
        'phi_take', CFG.vwa_phi_take, ...
        'compute_eta', false, ...
        'compute_phi_s', true, ...
        'small_kd_min', CFG.vwa_small_kd_min);

    zeros_field = zeros(size(eta11(:)));
    out = struct();
    out.order1 = zeros_field;
    out.order2 = vwa_out.phi_s{2}(:);
    out.order3 = vwa_out.phi_s{3}(:);
    out.order4 = zeros_field;
    out.meta = vwa_out.meta;
end

function spec_out = limit_linear_components_local(spec_in, max_components)
    spec_out = spec_in;
    if nargin < 2 || isempty(max_components) || ~isfinite(max_components)
        return;
    end

    max_components = max(1, round(max_components));
    if spec_in.num_components <= max_components
        return;
    end

    energy = abs(spec_in.coeff(:)).^2;
    [~, order] = sort(energy, 'descend');
    keep = sort(order(1:max_components));

    spec_out.indices = spec_in.indices(keep);
    spec_out.k = spec_in.k(keep);
    spec_out.coeff = spec_in.coeff(keep);
    spec_out.amplitude = spec_in.amplitude(keep);
    spec_out.phase = spec_in.phase(keep);
    spec_out.a = spec_in.a(keep);
    spec_out.b = spec_in.b(keep);
    spec_out.num_components = numel(keep);
end

function out = filter_surface_orders_x_only_local(in, x_vec, kp)
    out = in;
    out.eta1 = frequency_filtering_1d_local(in.eta1, x_vec, kp, 1);
    out.eta2 = frequency_filtering_1d_local(in.eta2, x_vec, kp, 2);
    out.eta3 = frequency_filtering_1d_local(in.eta3, x_vec, kp, 3);
    out.eta_nonlinear = out.eta1 + out.eta2 + out.eta3;
    order_names = {'order1', 'order2', 'order3'};
    for order_idx = 1:numel(order_names)
        field_name = order_names{order_idx};
        out.(field_name) = frequency_filtering_1d_local(in.(field_name), x_vec, kp, order_idx);
    end
    out.nonlinear = out.order1 + out.order2 + out.order3;
end

function out = build_surface_comparison_components_local(source, x_vec, kp, CFG, source_type)
    out = struct();
    out.first = extract_surface_component_local(source.order1, x_vec, kp, CFG, 'first');
    out.second_super = extract_surface_component_local(source.order2, x_vec, kp, CFG, 'second_super');
    out.third = extract_surface_component_local(source.order3, x_vec, kp, CFG, 'third');

    switch lower(source_type)
        case 'ow3d'
            out.second_sub = extract_surface_component_local(source.order4, x_vec, kp, CFG, 'second_sub');
        case 'mf12'
            out.second_sub = extract_surface_component_local(source.order2, x_vec, kp, CFG, 'second_sub');
        case 'vwa'
            out.second_sub = zeros(size(out.second_super));
        otherwise
            error('Unsupported source_type: %s', source_type);
    end
end

function field_out = extract_surface_component_local(field_in, x_vec, kp, CFG, component_name)
    field_in = field_in(:);
    switch lower(component_name)
        case 'first'
            band_factors = CFG.filter_band_order1;
        case 'second_super'
            band_factors = CFG.filter_band_order2;
        case 'third'
            band_factors = CFG.filter_band_order3;
        case 'second_sub'
            band_factors = CFG.filter_band_subharmonic;
        otherwise
            error('Unsupported component_name: %s', component_name);
    end
    field_out = apply_component_band_local(field_in, x_vec, kp, band_factors, CFG);
end

function field_out = apply_component_band_local(field_in, x_vec, kp, band_factors, CFG)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = smooth_bandpass_mask_local(abs(kx), band_factors(1) * kp, band_factors(2) * kp, ...
        CFG.filter_transition_low * kp, CFG.filter_transition_high * kp);
    field_out = real(ifft(fft(field_in) .* mask));
end

function split = compute_surface_frequency_split_local(ow3d_surface, mf12_surface, x_vec, kp, lowpass_cutoff, lowpass_transition)
    split = struct();
    split.meta = struct( ...
        'lowpass_cutoff', lowpass_cutoff, ...
        'lowpass_transition', lowpass_transition, ...
        'lowpass_cutoff_over_kp', lowpass_cutoff / kp, ...
        'lowpass_transition_over_kp', lowpass_transition / kp);

    component_names = {'first', 'second_super', 'third', 'second_sub'};
    for comp_idx = 1:numel(component_names)
        field_name = component_names{comp_idx};

        ow3d_raw = ow3d_surface.(field_name);
        mf12_raw = mf12_surface.(field_name);
        ow3d_band = ow3d_raw;
        mf12_band = mf12_raw;
        ow3d_low = lowpass_wavenumber_component_local(ow3d_raw, x_vec, lowpass_cutoff, lowpass_transition);
        mf12_low = lowpass_wavenumber_component_local(mf12_raw, x_vec, lowpass_cutoff, lowpass_transition);

        split.(field_name) = struct( ...
            'ow3d_raw', ow3d_raw(:), ...
            'mf12_raw', mf12_raw(:), ...
            'ow3d_band', ow3d_band(:), ...
            'mf12_band', mf12_band(:), ...
            'ow3d_lowpass', ow3d_low(:), ...
            'mf12_lowpass', mf12_low(:), ...
            'metrics_raw', compute_metrics_local(ow3d_raw, mf12_raw), ...
            'metrics_band', compute_metrics_local(ow3d_band, mf12_band), ...
            'metrics_lowpass', compute_metrics_local(ow3d_low, mf12_low));
    end
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

function [x_plot_shifted, eta11_shifted, field1_shifted, field2_shifted, field3_shifted, field4_shifted, x_center] = ...
        recenter_fields_for_plot_local(x_plot, eta11, eta_phases, field1, field2, field3, field4)
    x_plot = x_plot(:);
    eta11 = eta11(:);
    field1 = field1(:);
    field2 = field2(:);
    field3 = field3(:);
    field4 = field4(:);

    window_field = eta_phases(1, :).';
    if max(abs(window_field)) < 1e-12
        window_field = eta11;
    end

    shift_idx = compute_window_shift_local(window_field);
    x_plot_shifted = x_plot;
    eta11_shifted = circshift(eta11, shift_idx);
    field1_shifted = circshift(field1, shift_idx);
    field2_shifted = circshift(field2, shift_idx);
    field3_shifted = circshift(field3, shift_idx);
    field4_shifted = circshift(field4, shift_idx);

    recentered_window_field = circshift(window_field(:), shift_idx);
    recentered_envelope = abs(hilbert(recentered_window_field));
    [~, recentered_peak_idx] = max(recentered_envelope);
    x_center = x_plot_shifted(recentered_peak_idx);
end

function field_shifted = shift_with_window_local(field_in, eta_phases, eta11)
    field_in = field_in(:);
    window_field = eta_phases(1, :).';
    if max(abs(window_field)) < 1e-12
        window_field = eta11(:);
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

function y_limits = compute_pairwise_ylimits_local(fields_a, fields_b)
    n_fields = numel(fields_a);
    y_limits = zeros(n_fields, 2);
    for i = 1:n_fields
        values = [fields_a{i}(:); fields_b{i}(:)];
        y_abs_max = max(abs(values));
        if y_abs_max == 0
            y_abs_max = 1;
        end
        padding = 0.08 * y_abs_max;
        y_limits(i, :) = [-y_abs_max - padding, y_abs_max + padding];
    end
end

function y_limits = compute_threeway_ylimits_local(fields_a, fields_b, fields_c)
    n_fields = numel(fields_a);
    y_limits = zeros(n_fields, 2);
    for i = 1:n_fields
        values = [fields_a{i}(:); fields_b{i}(:); fields_c{i}(:)];
        y_abs_max = max(abs(values));
        if y_abs_max == 0
            y_abs_max = 1;
        end
        padding = 0.08 * y_abs_max;
        y_limits(i, :) = [-y_abs_max - padding, y_abs_max + padding];
    end
end

function draw_comparison_panel_local(ax, x_plot, ref_field, cmp_field, panel_title, y_limits, y_label, x_limits)
    plot(ax, x_plot, ref_field, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D'); hold(ax, 'on');
    plot(ax, x_plot, cmp_field, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12');
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 10);
end

function draw_threeway_comparison_panel_local(ax, x_plot, field_a, field_b, field_c, panel_title, y_limits, y_label, x_limits, label_a, label_b, label_c)
    plot(ax, x_plot, field_a, 'k-', 'LineWidth', 1.8, 'DisplayName', label_a); hold(ax, 'on');
    plot(ax, x_plot, field_b, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', label_b);
    plot(ax, x_plot, field_c, '-.', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'DisplayName', label_c);
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 10);
end

function [k_pos, amp_a, amp_b, amp_c, amp_d, amp_e, amp_f] = compute_three_pair_spectra_local( ...
        field_a, field_b, field_c, field_d, field_e, field_f, x_vec)
    field_a = field_a(:);
    field_b = field_b(:);
    field_c = field_c(:);
    field_d = field_d(:);
    field_e = field_e(:);
    field_f = field_f(:);
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
    fd = fft(field_d) / nx;
    fe = fft(field_e) / nx;
    ff = fft(field_f) / nx;

    amp_a = single_sided_amplitude_local(fa, pos_idx, nx);
    amp_b = single_sided_amplitude_local(fb, pos_idx, nx);
    amp_c = single_sided_amplitude_local(fc, pos_idx, nx);
    amp_d = single_sided_amplitude_local(fd, pos_idx, nx);
    amp_e = single_sided_amplitude_local(fe, pos_idx, nx);
    amp_f = single_sided_amplitude_local(ff, pos_idx, nx);
    k_pos = abs(kx(pos_idx));
end

function [k_pos, amp_ref, amp_cmp] = compute_spectrum_pair_local(field_ref, field_cmp, x_vec)
    field_ref = field_ref(:);
    field_cmp = field_cmp(:);
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

    f_ref = fft(field_ref) / nx;
    f_cmp = fft(field_cmp) / nx;
    amp_ref = single_sided_amplitude_local(f_ref, pos_idx, nx);
    amp_cmp = single_sided_amplitude_local(f_cmp, pos_idx, nx);
    k_pos = abs(kx(pos_idx));
end

function amp = single_sided_amplitude_local(fhat, pos_idx, nx)
    amp = 2 * abs(fhat(pos_idx));
    amp(1) = abs(fhat(1));
    if mod(nx, 2) == 0
        amp(end) = abs(fhat(pos_idx(end)));
    end
end

function plot_spectrum_pair_local(ax, k_plot, amp_ow3d, amp_mf12, panel_title)
    plot(ax, k_plot, amp_ow3d, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D'); hold(ax, 'on');
    plot(ax, k_plot, amp_mf12, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlabel(ax, '$|k| / k_p$', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, 'Amplitude', 'Interpreter', 'tex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 10);
end

function plot_spectrum_threeway_local(ax, k_plot, amp_ow3d, amp_mf12, amp_vwa, panel_title)
    plot(ax, k_plot, amp_ow3d, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D'); hold(ax, 'on');
    plot(ax, k_plot, amp_mf12, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12');
    plot(ax, k_plot, amp_vwa, '-.', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'DisplayName', 'VWA');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlabel(ax, '$|k| / k_p$', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, 'Amplitude', 'Interpreter', 'tex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 10);
end

function draw_difference_panel_local(ax, x_plot, diff_field, panel_title, y_limits, y_label, x_limits)
    plot(ax, x_plot, diff_field, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D - MF12'); hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    ylabel(ax, y_label, 'Interpreter', 'tex', 'FontSize', 13);
    xlabel(ax, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'FontSize', 10);
end

function field_out = lowpass_wavenumber_component_local(field_in, x_vec, k_cutoff, transition)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = smooth_bandpass_mask_local(abs(kx), 0.0, k_cutoff, 0.0, transition);
    field_out = real(ifft(fft(field_in) .* mask));
end

function print_frequency_split_metrics_local(surface_split)
    component_names = {'first', 'second_super', 'third', 'second_sub'};
    for comp_idx = 1:numel(component_names)
        field_name = component_names{comp_idx};
        raw = surface_split.(field_name).metrics_raw;
        band = surface_split.(field_name).metrics_band;
        low = surface_split.(field_name).metrics_lowpass;
        fprintf('%s raw:      corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
            field_name, raw.corr, raw.rmse, raw.peak_ratio);
        fprintf('%s band:     corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
            field_name, band.corr, band.rmse, band.peak_ratio);
        fprintf('%s lowpass:  corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', ...
            field_name, low.corr, low.rmse, low.peak_ratio);
    end
end
