% postprocess_ow3d_boundkinematics.m
% Reconstruct first-, second-, and third-order OW3D bound-kinematic
% components from four phase-shifted simulations using the same
% Hilbert/four-phase workflow as the surface postprocessor.
%
% Workflow:
% 1) read OW3D kinematics and compute time derivatives / pressure-like field
% 2) export the selected raw snapshot to MAT
% 3) run harmonic reconstruction and downstream postprocessing

clc;
clear;
close all;

CFG = struct();

% -------------------- User configuration --------------------
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check3');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_5.0_Akp_006_kd8.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.kinematics_file_id = 1; % Kinematics01.bin
CFG.phit_mode = 'uncorrected'; % 'uncorrected' -> Dt*phi, 'sigma_corrected' -> Dt*phi - w*sigma*etat
CFG.time_index = []; % [] -> use default near-final frame. Positive -> index from start. Negative -> index from end.
CFG.default_time_index_from_end = 160; % Used only when time_index = [].
CFG.lambda = 225;
CFG.gravity = 9.81;
CFG.kp_depth = 0.0279; % Use h = kd / kp with kd parsed from CFG.folder_pattern.
CFG.variables_to_process = {'u', 'w', 'phi'};
CFG.apply_x_filter = true;
CFG.sigma_mode = 'surface'; % 'surface', 'index', or 'value'
CFG.sigma_index = [];
CFG.sigma_value = 0.0;
CFG.save_raw_mat = false;
CFG.save_processed_mat = false;
CFG.vwa_surface_compare_variables = {'u', 'w'};
CFG.apply_vwa_eta11_filter = false;
CFG.vwa_small_kd_cutoff = 0.3;
CFG.plot_window_lambda = 5.0;
CFG.compare_mf12_subharmonic_surface = true;
CFG.mf12_linear_energy_keep = 0.9999;
CFG.mf12_subharmonic_cutoff_factor = 1.2;
CFG.mf12_subharmonic_transition_factor = 1.2;
CFG.export_standard_figures = true;
CFG.export_subharmonic_spectrum_figures = true;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');

CFG.vwa_required_surface_variables = resolve_required_vwa_surface_variables(CFG.vwa_surface_compare_variables);
CFG.process_variables = resolve_required_process_variables(CFG.variables_to_process, ...
    CFG.vwa_surface_compare_variables, CFG.vwa_required_surface_variables);

setup_vwa_surface_paths();

% -------------------- Load four OW3D kinematics snapshots ----------------
data_by_phase = cell(1, numel(CFG.phi_shifts_deg));
time_index_by_phase = zeros(1, numel(CFG.phi_shifts_deg));

for idx = 1:numel(CFG.phi_shifts_deg)
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
    if ~isfolder(case_folder)
        error('Missing phase folder: %s', case_folder);
    end

    kin_path = resolve_kinematics_path(case_folder, CFG.kinematics_file_id);
    if ~isfile(kin_path)
        error('Missing OW3D kinematics file: %s', kin_path);
    end

    phase_data = read_ow3d_kinematics_snapshot(kin_path, CFG.phit_mode);
    assert_phase_compatibility(phase_data, data_by_phase, idx);
    data_by_phase{idx} = phase_data;
    time_index_by_phase(idx) = resolve_time_index(CFG, CFG.time_index, phase_data.it, numel(phase_data.t));
    fprintf('Loaded %s\n', kin_path);
end

if any(time_index_by_phase ~= time_index_by_phase(1))
    error('Resolved time indices are inconsistent across the four phase files.');
end

selected_time_index = time_index_by_phase(1);
ref = data_by_phase{1};
x_vec = ref.x(:, 1);
y_vec = ref.y(1, :);
sigma_vec = ref.sigma(:);
t_vec = ref.t(:);
t_selected = t_vec(selected_time_index);
case_kd = extract_kd_from_case_pattern(CFG.folder_pattern);
depth_value = case_kd / CFG.kp_depth;

fprintf('Using kinematics time index %d of %d (t = %.6f s)\n', ...
    selected_time_index, numel(t_vec), t_selected);
fprintf('Using depth h = %.6f m from kd = %.4f and kp = %.4f 1/m\n', ...
    depth_value, case_kd, CFG.kp_depth);

% -------------------- Save raw snapshot before postprocessing ----------- 
if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

raw_meta = struct();
raw_meta.data_root = CFG.data_root;
raw_meta.folder_pattern = CFG.folder_pattern;
raw_meta.phi_shifts_deg = CFG.phi_shifts_deg;
raw_meta.kinematics_file_id = CFG.kinematics_file_id;
raw_meta.phit_mode = CFG.phit_mode;
raw_meta.time_index = selected_time_index;
raw_meta.time_value = t_selected;
raw_meta.lambda = CFG.lambda;
raw_meta.gravity = CFG.gravity;
raw_meta.kd = case_kd;
raw_meta.kp_depth = CFG.kp_depth;
raw_meta.depth_value = depth_value;
raw_meta.depth_source = 'h = kd / kp_depth parsed from CFG.folder_pattern';
raw_meta.sigma = sigma_vec;
raw_meta.x = x_vec;
raw_meta.y = y_vec;

if CFG.save_raw_mat
    raw_snapshot = build_raw_phase_snapshot(data_by_phase, CFG.phi_shifts_deg, selected_time_index);
    save(fullfile(CFG.output_dir, sprintf('OW3D_boundkinematics_raw_tidx_%04d.mat', selected_time_index)), ...
        'raw_snapshot', 'raw_meta', '-v7.3');
end

eta_phases = zeros(numel(CFG.phi_shifts_deg), size(ref.eta, 2));
vars_phases = initialize_variable_phase_storage(CFG.process_variables, numel(CFG.phi_shifts_deg), ref);

for idx = 1:numel(CFG.phi_shifts_deg)
    phase_data = data_by_phase{idx};
    eta_phases(idx, :) = squeeze(phase_data.eta(selected_time_index, :, 1));

    for v_idx = 1:numel(CFG.process_variables)
        var_name = CFG.process_variables{v_idx};
        var_field = squeeze(phase_data.(var_name)(selected_time_index, :, :, 1));
        vars_phases.(var_name)(idx, :, :) = var_field;
    end
end

% -------------------- Four-phase harmonic reconstruction -----------------
four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

eta_harmonics = reconstruct_harmonics_1d(eta_phases, four_phase_coef);
var_harmonics = struct();

for v_idx = 1:numel(CFG.process_variables)
    var_name = CFG.process_variables{v_idx};
    var_harmonics.(var_name) = reconstruct_harmonics_xz(vars_phases.(var_name), four_phase_coef);
end

% -------------------- Optional x-direction harmonic cleanup --------------
kp = 2 * pi / CFG.lambda;
if CFG.apply_x_filter
    eta_harmonics = filter_harmonics_x_only(eta_harmonics, x_vec, kp);

    for v_idx = 1:numel(CFG.process_variables)
        var_name = CFG.process_variables{v_idx};
        var_harmonics.(var_name) = filter_harmonics_x_only(var_harmonics.(var_name), x_vec, kp);
    end
end

% -------------------- VWA-like surface-kinematic approximation -----------
vwa_surface = struct();
eta11_surface = squeeze(eta_harmonics(1, :));
if CFG.apply_vwa_eta11_filter
    eta11_surface = frequency_filtering_1d_local(eta11_surface, x_vec, kp, 1);
end
sigma_idx = resolve_sigma_index(CFG, sigma_vec);
sigma_value = sigma_vec(sigma_idx);

for v_idx = 1:numel(CFG.vwa_surface_compare_variables)
    var_name = lower(CFG.vwa_surface_compare_variables{v_idx});
    vwa_surface.(var_name) = compute_surface_vwa_quantity_dispatch( ...
        var_name, eta11_surface, x_vec, depth_value, CFG.gravity, CFG.vwa_small_kd_cutoff, kp);
end
for v_idx = 1:numel(CFG.vwa_required_surface_variables)
    var_name = lower(CFG.vwa_required_surface_variables{v_idx});
    if isfield(vwa_surface, var_name)
        continue;
    end
    vwa_surface.(var_name) = compute_surface_vwa_quantity_dispatch( ...
        var_name, eta11_surface, x_vec, depth_value, CFG.gravity, CFG.vwa_small_kd_cutoff, kp);
end

surface_subharmonic_compare = struct();
if CFG.compare_mf12_subharmonic_surface
    if ~all(isfield(vars_phases, {'u', 'w'}))
        error('MF12 second-subharmonic surface comparison requires u and w in vars_phases.');
    end

    subharmonic_cutoff = CFG.mf12_subharmonic_cutoff_factor * kp;
    subharmonic_transition = CFG.mf12_subharmonic_transition_factor * subharmonic_cutoff;
    ow3d_u_surface_mean = squeeze(mean(vars_phases.u(:, sigma_idx, :), 1));
    ow3d_w_surface_mean = squeeze(mean(vars_phases.w(:, sigma_idx, :), 1));

    surface_subharmonic_compare.ow3d.u = lowpass_wavenumber_component_local( ...
        ow3d_u_surface_mean(:), x_vec, subharmonic_cutoff, subharmonic_transition);
    surface_subharmonic_compare.ow3d.w = lowpass_wavenumber_component_local( ...
        ow3d_w_surface_mean(:), x_vec, subharmonic_cutoff, subharmonic_transition);
    surface_subharmonic_compare.mf12 = compute_mf12_second_subharmonic_surface( ...
        eta11_surface(:), x_vec, depth_value, CFG.gravity, CFG.mf12_linear_energy_keep);
    surface_subharmonic_compare.meta = struct( ...
        'sigma_idx', sigma_idx, ...
        'sigma_value', sigma_value, ...
        'kp', kp, ...
        'cutoff', subharmonic_cutoff, ...
        'transition', subharmonic_transition, ...
        'energy_keep', CFG.mf12_linear_energy_keep, ...
        'linear_component_count', numel(surface_subharmonic_compare.mf12.linear_indices), ...
        'mf12_mode', 'difference_terms_only');

    fprintf('MF12 second subharmonic: kept %d linear components (energy keep = %.3f).\n', ...
        surface_subharmonic_compare.meta.linear_component_count, CFG.mf12_linear_energy_keep);
end

% -------------------- Save processed snapshot ----------------------------
meta = struct();
meta.data_root = CFG.data_root;
meta.folder_pattern = CFG.folder_pattern;
meta.phi_shifts_deg = CFG.phi_shifts_deg;
meta.kinematics_file_id = CFG.kinematics_file_id;
meta.phit_mode = CFG.phit_mode;
meta.time_index = selected_time_index;
meta.time_value = t_selected;
meta.lambda = CFG.lambda;
meta.kd = case_kd;
meta.kp_depth = CFG.kp_depth;
meta.depth_value = depth_value;
meta.depth_source = 'h = kd / kp_depth parsed from CFG.folder_pattern';
meta.sigma = sigma_vec;
meta.x = x_vec;
meta.y = y_vec;

if CFG.save_processed_mat
    save(fullfile(CFG.output_dir, sprintf('OW3D_boundkinematics_tidx_%04d.mat', selected_time_index)), ...
        'eta_harmonics', 'var_harmonics', 'meta', 'surface_subharmonic_compare', '-v7.3');
end

% -------------------- Visualization -------------------------------------
x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / CFG.lambda;
x_limits = resolve_plot_xlim(x_plot, eta_harmonics(1, :), CFG.plot_window_lambda);
line_colors = [0.10 0.10 0.10; 0.80 0.26 0.18; 0.12 0.39 0.71; 0.55 0.16 0.51];

if CFG.export_standard_figures
for v_idx = 1:numel(CFG.variables_to_process)
    var_name = CFG.variables_to_process{v_idx};

    if ismember(lower(var_name), lower(CFG.vwa_surface_compare_variables))
        continue;
    end

    harmonics = var_harmonics.(var_name);

    var1 = squeeze(harmonics(1, :, :)).';
    var2 = squeeze(harmonics(2, :, :)).';
    var3 = squeeze(harmonics(3, :, :)).';
    var_nl = var1 + var2 + var3;

    fields_to_plot = {var1(:, sigma_idx), var2(:, sigma_idx), var3(:, sigma_idx), var_nl(:, sigma_idx)};
    y_limits = compute_shared_ylimits_1d(fields_to_plot);

    fig = create_publishable_figure([120 100 1450 880]);
    tile = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('OW3D bound kinematics: %s at \\sigma = %.3f (t index = %d, t = %.4f s)', ...
        variable_display_name(var_name), sigma_value, selected_time_index, t_selected), ...
        'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    titles = { ...
        sprintf('(a) First-order %s', variable_display_name(var_name)), ...
        sprintf('(b) Second-order %s', variable_display_name(var_name)), ...
        sprintf('(c) Third-order %s', variable_display_name(var_name)), ...
        sprintf('(d) Nonlinear %s', variable_display_name(var_name))};

    for panel_idx = 1:4
        ax = nexttile(tile);
        draw_profile_panel(ax, x_plot, fields_to_plot{panel_idx}, line_colors(panel_idx, :), titles{panel_idx}, y_limits(panel_idx, :), variable_axis_label(var_name), x_limits);
    end

    annotation(fig, 'textbox', [0.13 0.01 0.82 0.04], ...
        'String', sprintf('Kinematics file %d, sigma index %d, sigma = %.4f, y = %.4f m.', ...
        CFG.kinematics_file_id, sigma_idx, sigma_value, y_vec(1)), ...
        'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontName', 'Times New Roman', 'FontSize', 11);

    exportgraphics(fig, fullfile(CFG.output_dir, ...
        sprintf('OW3D_boundkinematics_%s_sigma_%03d_tidx_%04d.png', var_name, sigma_idx, selected_time_index)), ...
        'Resolution', 300);
end

for v_idx = 1:numel(CFG.vwa_surface_compare_variables)
    var_name = lower(CFG.vwa_surface_compare_variables{v_idx});
    if ~isfield(var_harmonics, var_name) || ~isfield(vwa_surface, var_name)
        continue;
    end

    ow3d_harmonics = { ...
        squeeze(var_harmonics.(var_name)(1, sigma_idx, :)), ...
        squeeze(var_harmonics.(var_name)(2, sigma_idx, :)), ...
        squeeze(var_harmonics.(var_name)(3, sigma_idx, :))};
    vwa_harmonics = { ...
        vwa_surface.(var_name).order1(:), ...
        vwa_surface.(var_name).order2(:), ...
        vwa_surface.(var_name).order3(:)};
    y_limits = compute_pairwise_ylimits(ow3d_harmonics, vwa_harmonics);

    fig = create_publishable_figure([140 100 1450 920]);
    tile = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('Surface %s: OW3D vs VWA-like approximation (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
        quantity_display_name(var_name), sigma_value, selected_time_index, t_selected), ...
        'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    order_titles = {'(a) First harmonic', '(b) Second harmonic', '(c) Third harmonic'};
    for n = 1:3
        ax = nexttile(tile);
        draw_comparison_panel(ax, x_plot, ow3d_harmonics{n}, vwa_harmonics{n}, order_titles{n}, ...
            y_limits(n, :), quantity_axis_label(var_name), x_limits);
    end

    annotation(fig, 'textbox', [0.13 0.01 0.82 0.04], ...
        'String', build_surface_compare_footer(var_name, CFG.apply_vwa_eta11_filter, depth_value, y_vec(1), sigma_value), ...
        'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontName', 'Times New Roman', 'FontSize', 11);

    exportgraphics(fig, fullfile(CFG.output_dir, ...
        sprintf('OW3D_boundkinematics_surface_%s_vwa_compare_sigma_%03d_tidx_%04d.png', var_name, sigma_idx, selected_time_index)), ...
        'Resolution', 300);
end
end

if CFG.compare_mf12_subharmonic_surface && CFG.export_subharmonic_spectrum_figures && ~isempty(fieldnames(surface_subharmonic_compare))
    for quantity_name = {'u', 'w'}
        var_name = quantity_name{1};
        ref_field = surface_subharmonic_compare.ow3d.(var_name);
        mf12_field = surface_subharmonic_compare.mf12.(var_name);
        [k_plot, ow3d_amp, mf12_amp] = compute_one_sided_spectrum(ref_field, mf12_field, x_vec, kp);

        fig = create_publishable_figure([140 110 1450 760]);
        tile = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tile, sprintf('Surface %s second subharmonic spectrum: OW3D vs MF12 (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
            quantity_display_name(var_name), sigma_value, selected_time_index, t_selected), ...
            'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

        ax1 = nexttile(tile);
        plot(ax1, k_plot, ow3d_amp, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D'); hold(ax1, 'on');
        plot(ax1, k_plot, mf12_amp, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12');
        grid(ax1, 'on'); box(ax1, 'on');
        set(ax1, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
        xlabel(ax1, '$|k| / k_p$', 'Interpreter', 'latex', 'FontSize', 13);
        ylabel(ax1, 'Amplitude', 'Interpreter', 'tex', 'FontSize', 13);
        title(ax1, 'Linear scale', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
        xlim(ax1, [0, max(3, max(k_plot))]);
        legend(ax1, 'Location', 'best', 'FontSize', 10);

        ax2 = nexttile(tile);
        semilogy(ax2, k_plot, max(ow3d_amp, 1e-16), 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D'); hold(ax2, 'on');
        semilogy(ax2, k_plot, max(mf12_amp, 1e-16), '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12');
        grid(ax2, 'on'); box(ax2, 'on');
        set(ax2, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
        xlabel(ax2, '$|k| / k_p$', 'Interpreter', 'latex', 'FontSize', 13);
        ylabel(ax2, 'Amplitude (log)', 'Interpreter', 'tex', 'FontSize', 13);
        title(ax2, 'Semilog scale', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
        xlim(ax2, [0, max(3, max(k_plot))]);

        annotation(fig, 'textbox', [0.13 0.01 0.82 0.04], ...
            'String', build_surface_subharmonic_footer(surface_subharmonic_compare.meta, y_vec(1)), ...
            'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 11);

        exportgraphics(fig, fullfile(CFG.output_dir, ...
            sprintf('OW3D_boundkinematics_surface_%s_mf12_subharmonic_spectrum_sigma_%03d_tidx_%04d.png', ...
            var_name, sigma_idx, selected_time_index)), 'Resolution', 300);
    end
end

if CFG.compare_mf12_subharmonic_surface && CFG.export_standard_figures && ~isempty(fieldnames(surface_subharmonic_compare))
    for quantity_name = {'u', 'w'}
        var_name = quantity_name{1};
        ref_field = surface_subharmonic_compare.ow3d.(var_name);
        mf12_field = surface_subharmonic_compare.mf12.(var_name);
        metrics = compare_series_metrics(ref_field, mf12_field);
        y_limits_pair = compute_pairwise_ylimits({ref_field}, {mf12_field});
        y_limits_diff = compute_pairwise_ylimits({ref_field - mf12_field}, {zeros(size(ref_field))});

        fig = create_publishable_figure([140 110 1450 760]);
        tile = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tile, sprintf('Surface %s second subharmonic: OW3D vs MF12 (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
            quantity_display_name(var_name), sigma_value, selected_time_index, t_selected), ...
            'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

        ax1 = nexttile(tile);
        draw_comparison_panel(ax1, x_plot, ref_field, mf12_field, 'Subharmonic comparison', ...
            y_limits_pair(1, :), quantity_axis_label(var_name), x_limits, 'OW3D', 'MF12');
        text(ax1, 0.02, 0.92, sprintf('corr = %.3f, RMSE = %.3e, peak ratio = %.3f', ...
            metrics.corr, metrics.rmse, metrics.peak_ratio), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

        ax2 = nexttile(tile);
        draw_difference_panel(ax2, x_plot, ref_field - mf12_field, ...
            'Difference (OW3D - MF12)', y_limits_diff(1, :), quantity_axis_label(var_name), x_limits);

        annotation(fig, 'textbox', [0.13 0.01 0.82 0.04], ...
            'String', build_surface_subharmonic_footer(surface_subharmonic_compare.meta, y_vec(1)), ...
            'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
            'FontName', 'Times New Roman', 'FontSize', 11);

        exportgraphics(fig, fullfile(CFG.output_dir, ...
            sprintf('OW3D_boundkinematics_surface_%s_mf12_subharmonic_sigma_%03d_tidx_%04d.png', ...
            var_name, sigma_idx, selected_time_index)), 'Resolution', 300);
    end
end

disp('OW3D bound-kinematics postprocessing complete.');

% -------------------- Local helper functions -----------------------------
function setup_vwa_surface_paths()
    helper_dir = fullfile(fileparts(mfilename('fullpath')), 'test functions for VWA Opensource');
    mf12_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
    if ~isfolder(helper_dir)
        error('Missing VWA helper directory: %s', helper_dir);
    end
    if ~isfolder(mf12_dir)
        error('Missing MF12 helper directory: %s', mf12_dir);
    end
    addpath(helper_dir);
    addpath(mf12_dir);
end

function out = compute_surface_vwa_quantity_dispatch(quantity_name, eta11, x_vec, depth, gravity, small_kd_cutoff, kp)
    quantity_name = lower(quantity_name);

    switch quantity_name
        case {'u', 'w'}
            opts = struct('analytic_side', 'neg', 'small_kd_min', small_kd_cutoff);
            out = vwa_compute_surface_quantity(eta11, x_vec, depth, gravity, quantity_name, opts);
            out.kx = out.meta.kx;
            out.eta11 = eta11(:);
        otherwise
            out = approximate_surface_quantity_vwa_like( ...
                quantity_name, eta11, x_vec, depth, gravity, small_kd_cutoff, kp);
    end
end

function data = read_ow3d_kinematics_snapshot(kin_path, phit_mode)
    [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = ...
        read_kinematics_file_local(kin_path, phit_mode); %#ok<ASGLU>

    data = struct();
    data.it = it;
    data.eta = eta;
    data.etat_m = etat_m;
    data.etatt_m = etatt_m;
    data.phi = phi;
    data.phit = phit_m;
    data.phit_m = phit_m;
    data.p = p_m;
    data.pressure = p_m;
    data.ut = ut_m;
    data.u = u;
    data.v = v;
    data.w = w;
    data.uz = uz;
    data.vz = vz;
    data.wz = wz;
    data.x = x;
    data.y = y;
    data.h = h;
    data.sigma = sigma;
    data.t = t;
end

function [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = read_kinematics_file_local(file_path, phit_mode)
    nbits = 32;
    compute_derivatives = true;

    if nbits == 32
        int_nbit = 'int';
    elseif nbits == 64
        int_nbit = 'int64';
    else
        error('Illegal value for nbits: %d', nbits);
    end

    fid = fopen(file_path, 'r', 'ieee-le');
    if fid < 0
        error('Could not open kinematics file: %s', file_path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fread(fid, 1, int_nbit);
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
    fread(fid, 2, int_nbit);

    nx = floor((xend - xbeg) / xstride) + 1;
    ny = floor((yend - ybeg) / ystride) + 1;
    nt = floor((tend - tbeg) / tstride) + 1;

    tmp = zeros(nx * ny * max(nz, 5), 1);

    tmp(1:5 * nx * ny) = fread(fid, 5 * nx * ny, 'double');
    fread(fid, 2, int_nbit);

    x = zeros(nx, ny);
    y = zeros(nx, ny);
    h = zeros(nx, ny);
    x(:) = tmp(1:5:5 * nx * ny);
    y(:) = tmp(2:5:5 * nx * ny);
    h(:) = tmp(3:5:5 * nx * ny); %#ok<NASGU>

    for i = 1:nz
        sigma(i) = fread(fid, 1, 'double');
    end
    fread(fid, 2, int_nbit);

    eta = zeros(nt, nx, ny);
    etax = zeros(nt, nx, ny);
    etay = zeros(nt, nx, ny);
    phi = zeros(nt, nz, nx, ny);
    w = zeros(nt, nz, nx, ny);
    u = zeros(nt, nz, nx, ny);
    uz = zeros(nt, nz, nx, ny);
    v = zeros(nt, nz, nx, ny);
    vz = zeros(nt, nz, nx, ny);
    wz = zeros(nt, nz, nx, ny);
    t = (0:nt - 1) * dt * tstride;

    it = 0;
    for it_read = 1:nt - 1
        tmp_eta = fread(fid, nx * ny, 'double');
        if numel(tmp_eta) < nx * ny
            it = it_read - 1;
            break;
        end
        eta(it_read, :) = tmp_eta;
        fread(fid, 2, int_nbit);

        tmp_etax = fread(fid, nx * ny, 'double');
        if numel(tmp_etax) < nx * ny
            it = it_read - 1;
            break;
        end
        etax(it_read, :) = tmp_etax; %#ok<NASGU>
        fread(fid, 2, int_nbit);

        tmp_etay = fread(fid, nx * ny, 'double');
        if numel(tmp_etay) < nx * ny
            it = it_read - 1;
            break;
        end
        etay(it_read, :) = tmp_etay; %#ok<NASGU>
        fread(fid, 2, int_nbit);

        tmp_phi = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_phi) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        phi(it_read, :) = tmp_phi;
        fread(fid, 2, int_nbit);

        tmp_u = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_u) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        u(it_read, :) = tmp_u;
        fread(fid, 2, int_nbit);

        tmp_v = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_v) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        v(it_read, :) = tmp_v;
        fread(fid, 2, int_nbit);

        tmp_w = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_w) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        w(it_read, :) = tmp_w;
        fread(fid, 2, int_nbit);

        tmp_wz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_wz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        wz(it_read, :) = tmp_wz;
        fread(fid, 2, int_nbit);

        tmp_uz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_uz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        uz(it_read, :) = tmp_uz;
        fread(fid, 2, int_nbit);

        tmp_vz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_vz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        vz(it_read, :) = tmp_vz;
        fread(fid, 2, int_nbit);

        it = it_read;
    end

    if it <= 0
        error('No complete stored kinematics time step could be read from %s', file_path);
    end

    if it < nt
        eta = eta(1:it, :, :);
        phi = phi(1:it, :, :, :);
        u = u(1:it, :, :, :);
        v = v(1:it, :, :, :);
        w = w(1:it, :, :, :);
        uz = uz(1:it, :, :, :);
        vz = vz(1:it, :, :, :);
        wz = wz(1:it, :, :, :);
        t = t(1:it);
        nt = it;
    end

    if compute_derivatives
        alpha = 2;
        r = 2 * alpha + 1;
        c = build_stencil_even_local(alpha, 1);
        dt_matrix = spdiags(ones(nt, 1) * c(:, alpha + 1)', -alpha:alpha, nt, nt);
        for j = 1:alpha
            dt_matrix(j, :) = 0;
            dt_matrix(j, 1:r) = c(:, j)';
            dt_matrix(nt - j + 1, :) = 0;
            dt_matrix(nt - j + 1, nt - r + 1:nt) = c(:, r - j + 1)';
        end
        dt_matrix = dt_matrix / dt;

        etat_m = zeros(nt, size(eta, 2), size(eta, 3));
        etatt_m = zeros(nt, size(eta, 2), size(eta, 3));
        phit_m = zeros(size(phi));
        p_m = zeros(size(phi));
        ut_m = zeros(size(phi));

        for idy = 1:ny
            etat = zeros(nt, nx);
            etatt = zeros(nt, nx);
            phit = zeros(nt, nz, nx);
            p = zeros(nt, nz, nx);
            ut = zeros(nt, nz, nx);

            for ip = 1:nx
                eta_col = eta(:, ip, idy);
                etat(:, ip) = dt_matrix * eta_col;
                etatt(:, ip) = dt_matrix * etat(:, ip);

                for j = 1:nz
                    phi_col = phi(:, j, ip, idy);
                    w_col = w(:, j, ip, idy);
                    u_col = u(:, j, ip, idy);
                    uz_col = uz(:, j, ip, idy);

                    switch lower(phit_mode)
                        case 'uncorrected'
                            phit(:, j, ip) = dt_matrix * phi_col;
                        case 'sigma_corrected'
                            phit(:, j, ip) = dt_matrix * phi_col - w_col .* sigma(j) .* etat(:, ip);
                        otherwise
                            error('Unsupported phit_mode: %s', phit_mode);
                    end
                    p(:, j, ip) = -(phit(:, j, ip) + 0.5 * (u_col.^2 + v(:, j, ip, idy).^2 + w_col.^2));
                    ut(:, j, ip) = dt_matrix * u_col - uz_col .* sigma(j) .* etat(:, ip);
                end
            end

            etat_m(:, :, idy) = etat;
            etatt_m(:, :, idy) = etatt;
            phit_m(:, :, :, idy) = phit;
            p_m(:, :, :, idy) = p;
            ut_m(:, :, :, idy) = ut;
        end
    else
        etat_m = 0;
        etatt_m = 0;
        phit_m = 0;
        p_m = 0;
        ut_m = 0;
    end
end

function fx = build_stencil_even_local(alpha, der)
    rank = 2 * alpha + 1;
    fx = zeros(rank, rank);

    for ip = 1:alpha
        mat = zeros(rank, rank);
        row = 1;
        for m = -ip + 1:rank - ip
            for n = 1:rank
                mat(row, n) = m^(n - 1) / factorial(n - 1);
            end
            row = row + 1;
        end
        minv = inv(mat);
        fx(:, ip) = minv(der + 1, :).';
    end

    mat = zeros(rank, rank);
    row = 1;
    for m = -alpha:alpha
        for n = 1:rank
            mat(row, n) = m^(n - 1) / factorial(n - 1);
        end
        row = row + 1;
    end
    minv = inv(mat);
    fx(:, alpha + 1) = minv(der + 1, :).';

    if mod(der, 2) == 0
        for ip = 1:alpha
            fx(:, rank - ip + 1) = flipud(fx(:, ip));
        end
    else
        for ip = 1:alpha
            fx(:, rank - ip + 1) = -flipud(fx(:, ip));
        end
    end
end

function kin_path = resolve_kinematics_path(case_folder, file_id)
    if file_id < 10
        file_name = sprintf('Kinematics0%d.bin', file_id);
    else
        file_name = sprintf('Kinematics%d.bin', file_id);
    end

    kin_path = fullfile(case_folder, file_name);
end

function time_index = resolve_time_index(CFG, cfg_time_index, n_times_valid, n_times_total)
    if isempty(cfg_time_index)
        time_index = n_times_valid - CFG.default_time_index_from_end + 1;
    elseif isscalar(cfg_time_index) && cfg_time_index < 0
        time_index = n_times_valid + cfg_time_index + 1;
    else
        time_index = cfg_time_index;
    end
    time_index = round(time_index(1));
    time_index = min(max(1, time_index), min(n_times_valid, n_times_total));
end

function raw_snapshot = build_raw_phase_snapshot(data_by_phase, phi_shifts_deg, selected_time_index)
    raw_snapshot = struct();

    for idx = 1:numel(phi_shifts_deg)
        phase_field = matlab.lang.makeValidName(sprintf('phase_%03d', phi_shifts_deg(idx)));
        phase_data = data_by_phase{idx};

        snapshot = struct();
        snapshot.phi_shift_deg = phi_shifts_deg(idx);
        snapshot.it = phase_data.it;
        snapshot.t = phase_data.t(selected_time_index);
        snapshot.x = phase_data.x;
        snapshot.y = phase_data.y;
        snapshot.h = phase_data.h;
        snapshot.sigma = phase_data.sigma;
        snapshot.eta = squeeze(phase_data.eta(selected_time_index, :, :));
        snapshot.etat = squeeze(phase_data.etat_m(selected_time_index, :, :));
        snapshot.etatt = squeeze(phase_data.etatt_m(selected_time_index, :, :));
        snapshot.phi = squeeze(phase_data.phi(selected_time_index, :, :, :));
        snapshot.phit = squeeze(phase_data.phit(selected_time_index, :, :, :));
        snapshot.p = squeeze(phase_data.p(selected_time_index, :, :, :));
        snapshot.ut = squeeze(phase_data.ut(selected_time_index, :, :, :));
        snapshot.u = squeeze(phase_data.u(selected_time_index, :, :, :));
        snapshot.v = squeeze(phase_data.v(selected_time_index, :, :, :));
        snapshot.w = squeeze(phase_data.w(selected_time_index, :, :, :));
        snapshot.uz = squeeze(phase_data.uz(selected_time_index, :, :, :));
        snapshot.vz = squeeze(phase_data.vz(selected_time_index, :, :, :));
        snapshot.wz = squeeze(phase_data.wz(selected_time_index, :, :, :));

        raw_snapshot.(phase_field) = snapshot;
    end
end

function vars_phases = initialize_variable_phase_storage(var_names, n_phases, ref)
    vars_phases = struct();

    for v_idx = 1:numel(var_names)
        var_name = var_names{v_idx};
        var_size = size(ref.(var_name));
        nz = var_size(2);
        nx = var_size(3);
        vars_phases.(var_name) = zeros(n_phases, nz, nx);
    end
end

function harmonics = reconstruct_harmonics_1d(fields_by_phase, coef)
    analytic_part = hilbert(fields_by_phase.').';
    all_fields = cat(1, real(fields_by_phase), -imag(analytic_part));
    harmonics = zeros(4, size(fields_by_phase, 2));

    for n = 1:4
        harmonics(n, :) = coef(n, :) * all_fields;
    end
end

function harmonics = reconstruct_harmonics_xz(fields_by_phase, coef)
    analytic_part = hilbert_x_only_all(fields_by_phase);
    all_fields = cat(1, real(fields_by_phase), -imag(analytic_part));
    harmonics = zeros(4, size(fields_by_phase, 2), size(fields_by_phase, 3));

    for n = 1:4
        weights = reshape(coef(n, :), [8, 1, 1]);
        harmonics(n, :, :) = sum(all_fields .* weights, 1);
    end
end

function analytic_fields = hilbert_x_only_all(fields_by_phase)
    [n_phases, nz, nx] = size(fields_by_phase);
    analytic_fields = complex(zeros(n_phases, nz, nx));

    for phase_idx = 1:n_phases
        for z_idx = 1:nz
            analytic_fields(phase_idx, z_idx, :) = hilbert(squeeze(fields_by_phase(phase_idx, z_idx, :)));
        end
    end
end

function assert_phase_compatibility(phase_data, data_by_phase, idx)
    if idx == 1
        return;
    end

    ref = data_by_phase{1};

    if isempty(ref)
        return;
    end

    if ~isequal(size(phase_data.x), size(ref.x)) || ~isequal(size(phase_data.y), size(ref.y))
        error('Kinematics grid mismatch between phase 0 and phase index %d.', idx - 1);
    end

    if ~isequal(size(phase_data.h), size(ref.h)) || any(abs(phase_data.h(:) - ref.h(:)) > 1e-12)
        error('Bathymetry/depth mismatch between phase 0 and phase index %d.', idx - 1);
    end

    if ~isequal(size(phase_data.eta), size(ref.eta))
        error('Kinematics eta array size mismatch between phase 0 and phase index %d.', idx - 1);
    end

    if ~isequal(size(phase_data.u), size(ref.u))
        error('Kinematics velocity array size mismatch between phase 0 and phase index %d.', idx - 1);
    end

    if numel(phase_data.sigma) ~= numel(ref.sigma) || any(abs(phase_data.sigma(:) - ref.sigma(:)) > 1e-12)
        error('Sigma grid mismatch between phase 0 and phase index %d.', idx - 1);
    end

    if numel(phase_data.t) ~= numel(ref.t) || any(abs(phase_data.t(:) - ref.t(:)) > 1e-12)
        error('Stored kinematics time vector mismatch between phase 0 and phase index %d.', idx - 1);
    end
end

function harmonics_out = filter_harmonics_x_only(harmonics_in, x_vec, kp)
    harmonics_out = harmonics_in;

    if ndims(harmonics_in) == 2
        for n = 1:size(harmonics_in, 1)
            harmonics_out(n, :) = frequency_filtering_1d_local(squeeze(harmonics_in(n, :)), x_vec, kp, n);
        end
        return;
    end

    for n = 1:size(harmonics_in, 1)
        for z_idx = 1:size(harmonics_in, 2)
            harmonics_out(n, z_idx, :) = frequency_filtering_1d_local(squeeze(harmonics_in(n, z_idx, :)), x_vec, kp, n);
        end
    end
end

function field_out = frequency_filtering_1d_local(field_in, x_vec, kp, n)
    x_vec = x_vec(:);
    field_in = field_in(:);

    if numel(x_vec) ~= numel(field_in)
        error('Dimension mismatch between field and x vector.');
    end

    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;

    sigma_k = 0.5 * kp;
    k_target = n * kp;
    mask = exp(-((abs(kx) - k_target).^2) / (2 * sigma_k^2));

    field_fft = fft(field_in);
    field_out = ifft(field_fft .* mask);
    if isreal(field_in)
        field_out = real(field_out);
    end
end

function sigma_idx = resolve_sigma_index(CFG, sigma_vec)
    switch lower(CFG.sigma_mode)
        case 'surface'
            [~, sigma_idx] = max(sigma_vec);
        case 'index'
            sigma_idx = min(max(1, CFG.sigma_index), numel(sigma_vec));
        case 'value'
            [~, sigma_idx] = min(abs(sigma_vec - CFG.sigma_value));
        otherwise
            error('Unsupported sigma_mode: %s', CFG.sigma_mode);
    end
end

function fig = create_publishable_figure(fig_position)
    fig = figure('Color', 'w', 'Position', fig_position, 'Renderer', 'painters');
end

function draw_profile_panel(ax, x_plot, y_plot, line_color, panel_title, y_limits, y_label, x_limits)
    plot(ax, x_plot, y_plot, 'Color', line_color, 'LineWidth', 1.8, 'HandleVisibility', 'off');
    hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    ax.LineWidth = 1.0;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.XMinorGrid = 'off';
    ax.YMinorGrid = 'off';
    ax.GridAlpha = 0.14;
    ax.GridColor = [0 0 0];
    ax.Layer = 'top';
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    xlabel(ax, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
end

function y_limits = compute_shared_ylimits_1d(fields_to_plot)
    n_fields = numel(fields_to_plot);
    y_limits = zeros(n_fields, 2);

    for i = 1:n_fields
        values = fields_to_plot{i};
        y_abs_max = max(abs(values(:)));
        if y_abs_max == 0
            y_abs_max = 1;
        end
        padding = 0.08 * y_abs_max;
        y_limits(i, :) = [-y_abs_max - padding, y_abs_max + padding];
    end
end

function y_limits = compute_pairwise_ylimits(fields_a, fields_b)
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

function label = variable_display_name(var_name)
    switch lower(var_name)
        case 'u'
            label = 'horizontal velocity';
        case 'v'
            label = 'transverse velocity';
        case 'w'
            label = 'vertical velocity';
        case 'phi'
            label = 'velocity potential';
        case 'ut'
            label = 'horizontal acceleration';
        case 'phit'
            label = 'potential time derivative';
        case {'p', 'pressure'}
            label = 'pressure-like Bernoulli field';
        otherwise
            label = var_name;
    end
end

function out = approximate_surface_quantity_vwa_like(quantity_name, eta11, x_vec, depth, gravity, small_kd_cutoff, kp)
    eta11 = eta11(:);
    x_vec = x_vec(:);

    if numel(eta11) ~= numel(x_vec)
        error('VWA-like surface-u approximation requires eta11 and x_vec with matching length.');
    end

    if numel(x_vec) < 2
        error('Need at least two x-points for VWA-like surface-u approximation.');
    end

    dx = mean(diff(x_vec));
    if any(abs(diff(x_vec) - dx) > 1e-10 * max(1, abs(dx)))
        error('x_vec must be approximately uniformly spaced for the VWA-like surface-u approximation.');
    end

    nx = numel(x_vec);
    kx = vwa_kxgrid_local(nx, dx);

    if strcmpi(quantity_name, 'phit')
        out = approximate_surface_phit_bulk_like(eta11, x_vec, depth, gravity, small_kd_cutoff, kp, kx);
        return;
    end

    eta_analytic = hilbert(eta11);
    eta_hat = fft(eta_analytic);

    [coeff1, phase1] = surface_quantity_transfer_coeff(quantity_name, 1, abs(kx), depth, gravity, small_kd_cutoff);
    [coeff2, phase2] = surface_quantity_transfer_coeff(quantity_name, 2, abs(kx), depth, gravity, small_kd_cutoff);
    [coeff3, phase3] = surface_quantity_transfer_coeff(quantity_name, 3, abs(kx), depth, gravity, small_kd_cutoff);

    kappa1 = ifft(eta_hat .* coeff1);
    kappa2 = ifft(eta_hat .* coeff2);
    kappa3 = ifft(eta_hat .* coeff3);

    out = struct();
    out.order1 = vwa_apply_phase_operator(kappa1, phase1);
    out.order2 = vwa_apply_phase_operator(eta_analytic .* kappa2, phase2);
    out.order3 = vwa_apply_phase_operator((eta_analytic .^ 2) .* kappa3, phase3);
    out.kx = kx;
    out.eta11 = eta11;
end

function out = approximate_surface_phit_bulk_like(eta11, x_vec, depth, gravity, small_kd_cutoff, kp, kx)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    kx = kx(:);

    k_abs = abs(kx);
    kd = k_abs .* depth;
    kd_safe = max(kd, 1e-12);
    sigma = tanh(kd_safe);
    omega = sqrt(gravity .* k_abs .* sigma);
    zero_mask = (k_abs <= 1e-12);

    eta_analytic = hilbert(eta11);
    eta_hat = fft(eta_analytic);

    coeff_phit1 = -(omega.^2 ./ max(k_abs, 1e-12)) .* coth(kd_safe);
    coeff_phitz1 = -omega.^2;
    coeff_phitzz1 = -k_abs .* omega.^2 .* coth(kd_safe);

    coeff_eta2 = k_abs .* (3 - sigma.^2) ./ (8 * sigma.^3);
    coeff_phit2_direct = -omega.^2 .* (3 * cosh(2 * kd_safe) ./ (4 * sinh(kd_safe).^4));
    coeff_phitz2 = -3 .* k_abs .* omega.^2 .* cosh(kd_safe) ./ sinh(kd_safe).^3;
    coeff_phit3_direct = (3 / 64) .* k_abs .* omega.^2 .* ...
        ((-11 + 2 * cosh(2 * kd_safe)) .* cosh(3 * kd_safe) ./ sinh(kd_safe).^7);

    coeff_list = {coeff_phit1, coeff_phitz1, coeff_phitzz1, coeff_eta2, coeff_phit2_direct, coeff_phitz2, coeff_phit3_direct};
    for i = 1:numel(coeff_list)
        coeff = coeff_list{i};
        coeff(~isfinite(coeff)) = 0;
        coeff(zero_mask) = 0;
        coeff(kd < small_kd_cutoff) = 0;
        coeff_list{i} = coeff;
    end
    coeff_phit1 = coeff_list{1};
    coeff_phitz1 = coeff_list{2};
    coeff_phitzz1 = coeff_list{3};
    coeff_eta2 = coeff_list{4};
    coeff_phit2_direct = coeff_list{5};
    coeff_phitz2 = coeff_list{6};
    coeff_phit3_direct = coeff_list{7};

    phit1 = real(ifft(eta_hat .* coeff_phit1));
    phitz1 = real(ifft(eta_hat .* coeff_phitz1));
    phitzz1 = real(ifft(eta_hat .* coeff_phitzz1));

    eta2 = real(eta_analytic .* ifft(eta_hat .* coeff_eta2));
    phit2_direct = real(eta_analytic .* ifft(eta_hat .* coeff_phit2_direct));
    phitz2 = real(eta_analytic .* ifft(eta_hat .* coeff_phitz2));
    phit3_direct = real((eta_analytic .^ 2) .* ifft(eta_hat .* coeff_phit3_direct));

    phit2_surface_corr = frequency_filtering_1d_local(eta11 .* phitz1, x_vec, kp, 2);
    phit3_term_b = frequency_filtering_1d_local(eta11 .* phitz2, x_vec, kp, 3);
    phit3_term_c = frequency_filtering_1d_local(eta2 .* phitz1, x_vec, kp, 3);
    phit3_term_d = 0.5 * frequency_filtering_1d_local((eta11 .^ 2) .* phitzz1, x_vec, kp, 3);

    out = struct();
    out.order1 = phit1;
    out.order2 = phit2_direct + phit2_surface_corr;
    out.order3 = phit3_direct + phit3_term_b + phit3_term_c + phit3_term_d;
    out.kx = kx;
    out.eta11 = eta11;
    out.debug = struct( ...
        'eta2', eta2, ...
        'phitz1', phitz1, ...
        'phitzz1', phitzz1, ...
        'phit2_direct', phit2_direct, ...
        'phitz2', phitz2, ...
        'phit3_direct', phit3_direct, ...
        'phit2_surface_corr', phit2_surface_corr, ...
        'phit3_term_b', phit3_term_b, ...
        'phit3_term_c', phit3_term_c, ...
        'phit3_term_d', phit3_term_d);
end

function out = compute_mf12_second_subharmonic_surface(eta11, x_vec, depth, gravity, energy_keep)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    nx = numel(x_vec);

    if nx < 2
        error('MF12 second-order surface calculation requires at least two x-points.');
    end

    dx = x_vec(2) - x_vec(1);
    kx_grid = vwa_kxgrid_local(nx, dx);
    eta_hat = fft(eta11) / nx;

    if mod(nx, 2) == 0
        positive_idx = 2:(nx / 2);
    else
        positive_idx = 2:((nx + 1) / 2);
    end

    positive_idx = positive_idx(abs(eta_hat(positive_idx)) > 1e-12 * max(abs(eta_hat)));
    positive_idx = select_energy_dominant_indices_local(eta_hat, positive_idx, energy_keep);

    out = struct();
    if isempty(positive_idx)
        out.u = zeros(nx, 1);
        out.w = zeros(nx, 1);
        out.phi = zeros(nx, 1);
        out.eta = zeros(nx, 1);
        out.linear_indices = [];
        out.energy_keep = energy_keep;
        return;
    end

    kx = kx_grid(positive_idx).';
    ky = zeros(size(kx));
    a = 2 * real(eta_hat(positive_idx)).';
    b = 2 * imag(eta_hat(positive_idx)).';

    coeffs2 = mf12_direct_coefficients(2, gravity, depth, a, b, kx, ky, 0, 0, 0);
    [eta20, phi20, u20, ~, w20] = mf12_second_subharmonic_kinematics(coeffs2, x_vec.', 0, 0, 0);

    out.eta = eta20(:);
    out.u = u20(:);
    out.w = w20(:);
    out.phi = phi20(:);
    out.linear_indices = positive_idx(:);
    out.energy_keep = energy_keep;
end

function field_out = lowpass_wavenumber_component_local(field_in, x_vec, k_cutoff, transition)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = exp(-(abs(kx) / max(transition, dkx)).^4);
    mask(abs(kx) <= k_cutoff) = 1;
    field_out = real(ifft(fft(field_in) .* mask));
end

function keep_idx = select_energy_dominant_indices_local(spectrum, candidate_idx, energy_keep)
    if nargin < 3 || isempty(energy_keep)
        energy_keep = 1.0;
    end

    energy_keep = min(max(energy_keep, 0), 1);
    if isempty(candidate_idx)
        keep_idx = candidate_idx;
        return;
    end

    spectral_energy = abs(spectrum(candidate_idx)).^2;
    total_energy = sum(spectral_energy);
    if total_energy <= 0 || energy_keep >= 1
        keep_idx = candidate_idx;
        return;
    end

    [sorted_energy, order] = sort(spectral_energy, 'descend');
    cumulative_energy = cumsum(sorted_energy) / total_energy;
    cutoff_idx = find(cumulative_energy >= energy_keep, 1, 'first');
    keep_idx = sort(candidate_idx(order(1:cutoff_idx)));
end

function [coeff, phase_type] = surface_quantity_transfer_coeff(quantity_name, order, k_abs, depth, gravity, small_kd_cutoff)
    if ismember(lower(quantity_name), {'u', 'w'})
        [coeff, phase_type] = vwa_surface_quantity_coeff(quantity_name, order, k_abs, depth, gravity, small_kd_cutoff);
        return;
    end

    kd = k_abs .* depth;
    kd_safe = max(kd, 1e-12);
    sigma = tanh(kd_safe);
    omega = sqrt(gravity .* k_abs .* sigma);

    coeff = zeros(size(k_abs));
    zero_mask = (k_abs <= 1e-12);
    phase_type = 'real';

    switch lower(quantity_name)
        case 'u'
            phase_type = 'real';
            switch order
                case 1
                    coeff(~zero_mask) = omega(~zero_mask) .* coth(kd_safe(~zero_mask));
                case 2
                    coeff(~zero_mask) = k_abs(~zero_mask) .* omega(~zero_mask) .* ...
                        (0.5 + 3 * cosh(2 * kd_safe(~zero_mask)) ./ (4 * sinh(kd_safe(~zero_mask)).^4));
                case 3
                    term1 = -(3 / 64) * k_abs.^2 .* omega .* ...
                        ((-11 + 2 * cosh(2 * kd_safe)) ./ sinh(kd_safe).^7) .* cosh(3 * kd_safe);
                    term2 = (3 / 2) * k_abs.^2 .* omega .* cosh(kd_safe) ./ sinh(kd_safe).^3;
                    term3 = (1 / 8) * k_abs.^2 .* omega .* ...
                        (2 + cosh(2 * kd_safe)) .* cosh(kd_safe) ./ sinh(kd_safe).^3;
                    term4 = (1 / 8) * omega .* k_abs.^2 .* coth(kd_safe);
                    coeff(~zero_mask) = term1(~zero_mask) + term2(~zero_mask) + term3(~zero_mask) + term4(~zero_mask);
                otherwise
                    error('Unsupported order %d for quantity %s.', order, quantity_name);
            end

        case 'w'
            phase_type = 'neg_imag';
            switch order
                case 1
                    coeff(~zero_mask) = -omega(~zero_mask);
                case 2
                    coeff(~zero_mask) = -omega(~zero_mask) .* k_abs(~zero_mask) .* ...
                        (0.5 .* coth(kd_safe(~zero_mask)) + 1.5 .* cosh(kd_safe(~zero_mask)) ./ sinh(kd_safe(~zero_mask)).^3);
                case 3
                    term1 = (3 / 64) * omega .* k_abs.^2 .* ...
                        ((-11 + 2 * cosh(2 * kd_safe)) ./ sinh(kd_safe).^7) .* sinh(3 * kd_safe);
                    term2 = -(3 / 4) * omega .* k_abs.^2 .* cosh(2 * kd_safe) ./ sinh(kd_safe).^4;
                    term3 = -(1 / 8) * omega .* k_abs.^2 .* ...
                        (2 + cosh(2 * kd_safe)) .* coth(kd_safe).^2 ./ sinh(kd_safe).^2;
                    term4 = -(1 / 8) * omega .* k_abs.^2;
                    coeff(~zero_mask) = term1(~zero_mask) + term2(~zero_mask) + term3(~zero_mask) + term4(~zero_mask);
                otherwise
                    error('Unsupported order %d for quantity %s.', order, quantity_name);
            end

        case 'phi'
            phase_type = 'imag';
            switch order
                case 1
                    coeff(~zero_mask) = -(omega(~zero_mask) ./ k_abs(~zero_mask)) .* coth(kd_safe(~zero_mask));
                case 2
                    coeff(~zero_mask) = -(omega(~zero_mask) / 8) .* ...
                        (4 + 3 * coth(kd_safe(~zero_mask)) ./ sinh(kd_safe(~zero_mask)).^3);
                case 3
                    coeff(~zero_mask) = -(k_abs(~zero_mask) ./ (64 * omega(~zero_mask))) .* ...
                        (8 * gravity * k_abs(~zero_mask) + omega(~zero_mask).^2 .* coth(kd_safe(~zero_mask)) .* ...
                        (16 + 56 ./ sinh(kd_safe(~zero_mask)).^2 + 32 ./ sinh(kd_safe(~zero_mask)).^4 + ...
                        9 ./ sinh(kd_safe(~zero_mask)).^6));
                otherwise
                    error('Unsupported order %d for quantity %s.', order, quantity_name);
            end

        case 'phit'
            phase_type = 'real';
            switch order
                case 1
                    coeff(~zero_mask) = -(omega(~zero_mask).^2 ./ k_abs(~zero_mask)) .* coth(kd_safe(~zero_mask));
                case 2
                    coeff(~zero_mask) = -omega(~zero_mask).^2 .* ...
                        (0.5 + 3 * cosh(2 * kd_safe(~zero_mask)) ./ (4 * sinh(kd_safe(~zero_mask)).^4));
                case 3
                    term1 = (3 / 64) .* ((-11 + 2 * cosh(2 * kd_safe)) .* cosh(3 * kd_safe) ./ sinh(kd_safe).^7);
                    term2 = -(3 / 2) .* cosh(kd_safe) ./ sinh(kd_safe).^3;
                    term3 = -((2 + cosh(2 * kd_safe)) .* coth(kd_safe)) ./ (8 * sinh(kd_safe).^2);
                    term4 = -(1 / 8) .* coth(kd_safe);
                    coeff(~zero_mask) = k_abs(~zero_mask) .* omega(~zero_mask).^2 .* ...
                        (term1(~zero_mask) + term2(~zero_mask) + term3(~zero_mask) + term4(~zero_mask));
                otherwise
                    error('Unsupported order %d for quantity %s.', order, quantity_name);
            end

        otherwise
            error('Unsupported surface quantity for VWA-like approximation: %s', quantity_name);
    end

    coeff(~isfinite(coeff)) = 0;
    coeff(zero_mask) = 0;
    coeff(kd < small_kd_cutoff) = 0;
end

function required_vars = resolve_required_vwa_surface_variables(compare_vars)
    compare_vars = lower(compare_vars(:).');
    required_vars = compare_vars;
end

function process_vars = resolve_required_process_variables(base_vars, compare_vars, required_surface_vars)
    base_vars = lower(base_vars(:).');
    compare_vars = lower(compare_vars(:).');
    required_surface_vars = lower(required_surface_vars(:).');
    process_vars = unique([base_vars, compare_vars, required_surface_vars], 'stable');
end

function kx = vwa_kxgrid_local(nx, dx)
    if nx < 2
        kx = 0;
        return;
    end

    dk = 2 * pi / (nx * dx);
    kpos = 0:floor(nx / 2);
    kneg = -ceil(nx / 2) + 1:-1;
    kx = (dk * [kpos, kneg]).';
end

function draw_comparison_panel(ax, x_plot, y_ow3d, y_vwa, panel_title, y_limits, y_label, x_limits, label_a, label_b)
    if nargin < 9 || isempty(label_a)
        label_a = 'OW3D';
    end
    if nargin < 10 || isempty(label_b)
        label_b = 'VWA-like';
    end

    h1 = plot(ax, x_plot, y_ow3d, 'Color', [0.10 0.10 0.10], 'LineWidth', 1.8, 'DisplayName', 'OW3D');
    hold(ax, 'on');
    h2 = plot(ax, x_plot, y_vwa, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', label_b);
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    ax.LineWidth = 1.0;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.GridAlpha = 0.14;
    ax.GridColor = [0 0 0];
    ax.Layer = 'top';
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    xlabel(ax, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
    h1.DisplayName = label_a;
    h2.DisplayName = label_b;
    legend(ax, [h1, h2], {label_a, label_b}, 'Location', 'best', 'FontSize', 10);
end

function draw_difference_panel(ax, x_plot, y_diff, panel_title, y_limits, y_label, x_limits)
    plot(ax, x_plot, y_diff, 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'HandleVisibility', 'off');
    hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9, 'HandleVisibility', 'off');
    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    ax.LineWidth = 1.0;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.GridAlpha = 0.14;
    ax.GridColor = [0 0 0];
    ax.Layer = 'top';
    xlim(ax, x_limits);
    ylim(ax, y_limits);
    xlabel(ax, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, y_label, 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
end

function metrics = compare_series_metrics(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    metrics = struct('corr', NaN, 'rmse', NaN, 'peak_ratio', NaN);

    if isempty(reference) || isempty(candidate)
        return;
    end

    if ~(all(abs(reference) < eps) || all(abs(candidate) < eps))
        cc = corrcoef(reference, candidate);
        metrics.corr = cc(1, 2);
    end
    metrics.rmse = sqrt(mean((reference - candidate).^2));

    ref_peak = max(abs(reference));
    if ref_peak > eps
        metrics.peak_ratio = max(abs(candidate)) / ref_peak;
    end
end

function [k_plot, amp_a, amp_b] = compute_one_sided_spectrum(field_a, field_b, x_vec, kp)
    field_a = field_a(:);
    field_b = field_b(:);
    x_vec = x_vec(:);
    nx = numel(x_vec);
    dx = x_vec(2) - x_vec(1);
    kx = vwa_kxgrid_local(nx, dx);

    fft_a = fft(field_a) / nx;
    fft_b = fft(field_b) / nx;

    if mod(nx, 2) == 0
        positive_idx = 1:(nx / 2 + 1);
    else
        positive_idx = 1:((nx + 1) / 2);
    end

    k_plot = abs(kx(positive_idx)) / kp;
    amp_a = abs(fft_a(positive_idx));
    amp_b = abs(fft_b(positive_idx));

    if numel(positive_idx) > 2
        amp_a(2:end-1) = 2 * amp_a(2:end-1);
        amp_b(2:end-1) = 2 * amp_b(2:end-1);
    end
end

function footer = build_surface_subharmonic_footer(meta, y_value)
    footer = sprintf(['OW3D reference: four-phase mean + low-pass. MF12 input: reconstructed $\\eta^{(1)}$. ' ...
        '$|k|<%.2f k_p$, transition = %.2f k_p, retained linear components = %d, y = %.4f m.'], ...
        meta.cutoff / meta.kp, meta.transition / meta.kp, meta.linear_component_count, y_value);
end

function kd = extract_kd_from_case_pattern(folder_pattern)
    token = regexp(folder_pattern, 'kd(?<kd>\d+(?:\.\d+)?)', 'names', 'once');
    if isempty(token) || ~isfield(token, 'kd')
        error('Unable to parse kd from CFG.folder_pattern: %s', folder_pattern);
    end

    kd = str2double(token.kd);
    if ~(isfinite(kd) && kd > 0)
        error('Parsed invalid kd value from CFG.folder_pattern: %s', folder_pattern);
    end
end

function x_limits = resolve_plot_xlim(x_plot, eta11, plot_window_lambda)
    eta11 = eta11(:);
    x_plot = x_plot(:);

    if numel(eta11) ~= numel(x_plot)
        x_limits = [x_plot(1), x_plot(end)];
        return;
    end

    [~, peak_idx] = max(abs(eta11));
    half_width = 0.5 * plot_window_lambda;
    x_center = x_plot(peak_idx);
    x_limits = [x_center - half_width, x_center + half_width];
    x_limits(1) = max(x_limits(1), x_plot(1));
    x_limits(2) = min(x_limits(2), x_plot(end));
end

function out = ternary_text(condition, true_text, false_text)
    if condition
        out = true_text;
    else
        out = false_text;
    end
end

function footer_text = build_surface_compare_footer(var_name, is_eta_filtered, depth_value, y_value, sigma_value)
    quantity_note = sprintf('VWA-like input uses %s OW3D linear surface elevation \\eta^{(1)} only.', ...
        ternary_text(is_eta_filtered, 'filtered', 'unfiltered'));
    footer_text = sprintf('%s Depth = %.4f m, y = %.4f m, \\sigma = %.4f.', ...
        quantity_note, depth_value, y_value, sigma_value);
end

function label = quantity_display_name(var_name)
    switch lower(var_name)
        case 'u'
            label = 'horizontal velocity';
        case 'w'
            label = 'vertical velocity';
        case 'phi'
            label = 'surface potential';
        case 'phit'
            label = 'surface potential time derivative';
        case {'p', 'pressure'}
            label = 'surface pressure-like field';
        otherwise
            label = var_name;
    end
end

function label = quantity_axis_label(var_name)
    switch lower(var_name)
        case {'u', 'w'}
            label = sprintf('$%s_s$ (m/s)', var_name);
        case 'phi'
            label = '$\phi_s$ (m$^2$/s)';
        case 'phit'
            label = '$\phi_{s,t}$ (m$^2$/s$^2$)';
        case {'p', 'pressure'}
            label = '$p_s$ (m$^2$/s$^2$)';
        otherwise
            label = ['$', var_name, '$'];
    end
end

function label = variable_axis_label(var_name)
    switch lower(var_name)
        case {'u', 'v', 'w', 'ut'}
            label = sprintf('$%s$', var_name);
        case 'phi'
            label = '$\phi$';
        case 'phit'
            label = '$\phi_t$';
        case {'p', 'pressure'}
            label = '$p$';
        otherwise
            label = ['$', var_name, '$'];
    end
end
