% compare_ow3d_mf12_subharmonic_amplitude_budget.m
% Unified amplitude-budget diagnostics for OW3D vs MF12 second-order
% subharmonic u20 and w20.

clc;
clear;
close all;

CFG = struct();
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check3');
CFG.case_patterns = { ...
    'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d', ...
    'T_init-20_Tp_Alpha_5.0_Akp_006_kd8.0_phi_%d'};
CFG.phi_shifts_deg = 0:90:270;
CFG.kinematics_file_id = 1;
CFG.phit_mode = 'uncorrected';
CFG.time_index = [];
CFG.default_time_index_from_end = 160;
CFG.lambda = 225;
CFG.gravity = 9.81;
CFG.kp_depth = 0.0279;
CFG.apply_x_filter = true;
CFG.sigma_mode = 'surface';
CFG.sigma_index = [];
CFG.sigma_value = 0.0;
CFG.linear_energy_keep = 0.99999;
CFG.local_window_lambda = 6.0;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');

folder_pattern_env = getenv('OW3D_FOLDER_PATTERN');
if ~isempty(folder_pattern_env)
    CFG.case_patterns = {folder_pattern_env};
end

mf12_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
if ~isfolder(mf12_dir)
    error('Missing MF12 source directory: %s', mf12_dir);
end
addpath(mf12_dir);

for case_idx = 1:numel(CFG.case_patterns)
    case_pattern = CFG.case_patterns{case_idx};
    process_case_local(CFG, case_pattern);
end

disp('Subharmonic amplitude-budget diagnostics complete.');

function process_case_local(CFG, folder_pattern)
    data_by_phase = cell(1, numel(CFG.phi_shifts_deg));
    time_index_by_phase = zeros(1, numel(CFG.phi_shifts_deg));

    for idx = 1:numel(CFG.phi_shifts_deg)
        case_folder = fullfile(CFG.data_root, sprintf(folder_pattern, CFG.phi_shifts_deg(idx)));
        kin_path = resolve_kinematics_path_local(case_folder, CFG.kinematics_file_id);
        phase_data = read_ow3d_kinematics_snapshot_local(kin_path, CFG.phit_mode);
        assert_phase_compatibility_local(phase_data, data_by_phase, idx);
        data_by_phase{idx} = phase_data;
        time_index_by_phase(idx) = resolve_time_index_local(CFG, CFG.time_index, phase_data.it, numel(phase_data.t));
        fprintf('Loaded %s\n', kin_path);
    end

    if any(time_index_by_phase ~= time_index_by_phase(1))
        error('Resolved time indices are inconsistent across the four phase files.');
    end

    selected_time_index = time_index_by_phase(1);
    ref = data_by_phase{1};
    x_vec = ref.x(:, 1);
    sigma_vec = ref.sigma(:);
    t_vec = ref.t(:);
    t_selected = t_vec(selected_time_index);
    case_kd = extract_kd_from_case_pattern_local(folder_pattern);
    depth_value = case_kd / CFG.kp_depth;
    kp = 2 * pi / CFG.lambda;
    case_folder_ref = fullfile(CFG.data_root, sprintf(folder_pattern, CFG.phi_shifts_deg(1)));
    ow3d_diff_alpha = resolve_ow3d_diff_alpha_local(case_folder_ref, []);
    case_tag = build_case_tag_local(folder_pattern);
    output_dir = fullfile(CFG.output_dir, case_tag);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    fprintf('Using kinematics time index %d of %d (t = %.6f s)\n', ...
        selected_time_index, numel(t_vec), t_selected);
    fprintf('Using depth h = %.6f m from kd = %.4f and kp = %.4f 1/m\n', ...
        depth_value, case_kd, CFG.kp_depth);
    fprintf('Using OW3D DiffXEven alpha = %d\n', ow3d_diff_alpha);

    eta_phases = zeros(numel(CFG.phi_shifts_deg), size(ref.eta, 2));
    vars_phases = struct();
    vars_phases.u = zeros(numel(CFG.phi_shifts_deg), size(ref.u, 2), size(ref.u, 3));
    vars_phases.w = zeros(numel(CFG.phi_shifts_deg), size(ref.w, 2), size(ref.w, 3));
    vars_phases.phi = zeros(numel(CFG.phi_shifts_deg), size(ref.phi, 2), size(ref.phi, 3));

    for idx = 1:numel(CFG.phi_shifts_deg)
        phase_data = data_by_phase{idx};
        eta_phases(idx, :) = squeeze(phase_data.eta(selected_time_index, :, 1));
        vars_phases.u(idx, :, :) = squeeze(phase_data.u(selected_time_index, :, :, 1));
        vars_phases.w(idx, :, :) = squeeze(phase_data.w(selected_time_index, :, :, 1));
        vars_phases.phi(idx, :, :) = squeeze(phase_data.phi(selected_time_index, :, :, 1));
    end

    four_phase_coef = [
        0.25  0    -0.25  0     0    -0.25  0     0.25;
        0.25 -0.25  0.25 -0.25  0     0      0     0;
        0.25  0    -0.25  0     0     0.25  0    -0.25;
        0.25  0.25  0.25  0.25  0     0      0     0];

    eta_harmonics = reconstruct_harmonics_1d_local(eta_phases, four_phase_coef);
    if CFG.apply_x_filter
        eta_harmonics = filter_harmonics_x_only_local(eta_harmonics, x_vec, kp);
    end

    sigma_idx = resolve_sigma_index_local(CFG, sigma_vec);
    sigma_value = sigma_vec(sigma_idx);
    eta11_surface = squeeze(eta_harmonics(1, :)).';

    ow3d = compute_ow3d_subharmonic_targets_local( ...
        eta_phases, vars_phases, ref.h(:, 1), sigma_idx, sigma_value, x_vec, four_phase_coef, ow3d_diff_alpha);
    mf12 = compute_mf12_subharmonic_targets_local(eta11_surface, x_vec, depth_value, CFG.gravity, CFG.linear_energy_keep);

    window_info = select_local_window_local(eta11_surface, x_vec, CFG.lambda, CFG.local_window_lambda, sigma_idx);

    metrics = struct();
    metrics.u20 = struct();
    metrics.u20.mf12_u20 = compare_with_gain_local(ow3d.u_bare, mf12.u20, window_info.mask);
    metrics.u20.bulk_phi_x = compare_with_gain_local(ow3d.u_bare, mf12.phix_bulk, window_info.mask);
    metrics.u20.surface_phi_x = compare_with_gain_local(ow3d.u_bare, mf12.phix_surface, window_info.mask);

    metrics.w20 = struct();
    metrics.w20.mf12_w20 = compare_with_gain_local(ow3d.w20, mf12.w20_z0, window_info.mask);
    metrics.w20.mf12_w20_z0 = compare_with_gain_local(ow3d.w20, mf12.w20_z0, window_info.mask);
    metrics.w20.mf12_w20_zeta11 = compare_with_gain_local(ow3d.w20, mf12.w20_zeta11, window_info.mask);
    metrics.w20.main_band_peak_ratio_z0 = compute_main_band_peak_ratio_local(ow3d.w20, mf12.w20_z0, x_vec, kp, 4);
    metrics.w20.main_band_peak_ratio_zeta11 = compute_main_band_peak_ratio_local(ow3d.w20, mf12.w20_zeta11, x_vec, kp, 4);

    fprintf('\n=== Subharmonic amplitude budget: %s ===\n', case_tag);
    print_metric_line_local('u20 : OW3D bare vs MF12 u20', metrics.u20.mf12_u20);
    print_metric_line_local('u20 : OW3D bare vs MF12 bulk phi20_x', metrics.u20.bulk_phi_x);
    print_metric_line_local('u20 : OW3D bare vs MF12 surface phi20_x', metrics.u20.surface_phi_x);
    print_metric_line_local('w20 : OW3D vs MF12 w20(z=0)', metrics.w20.mf12_w20_z0);
    print_metric_line_local('w20 : OW3D vs MF12 w20(z=eta11)', metrics.w20.mf12_w20_zeta11);
    fprintf('w20 main-band peak ratio (z=0)    : %.6f\n', metrics.w20.main_band_peak_ratio_z0);
    fprintf('w20 main-band peak ratio (z=eta11): %.6f\n', metrics.w20.main_band_peak_ratio_zeta11);

    create_amplitude_budget_figure_local(output_dir, case_tag, x_vec, CFG.lambda, sigma_value, ...
        selected_time_index, t_selected, ow3d, mf12, metrics, window_info);

    results = struct();
    results.case_tag = case_tag;
    results.t_index = selected_time_index;
    results.t_value = t_selected;
    results.sigma_idx = sigma_idx;
    results.sigma_value = sigma_value;
    results.window_info = window_info;
    results.ow3d = ow3d;
    results.mf12 = mf12;
    results.metrics = metrics;
    results.meta = struct( ...
        'folder_pattern', folder_pattern, ...
        'gravity', CFG.gravity, ...
        'lambda', CFG.lambda, ...
        'kp', kp, ...
        'depth', depth_value, ...
        'linear_energy_keep', CFG.linear_energy_keep, ...
        'ow3d_diff_alpha', ow3d_diff_alpha, ...
        'local_window_lambda', CFG.local_window_lambda);

    save(fullfile(output_dir, sprintf('compare_ow3d_mf12_subharmonic_amplitude_budget_sigma_%03d_tidx_%04d.mat', ...
        sigma_idx, selected_time_index)), 'results', '-v7.3');
end

function create_amplitude_budget_figure_local(output_dir, case_tag, x_vec, lambda, sigma_value, t_index, t_value, ow3d, mf12, metrics, window_info)
    x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / lambda;
    x_limits = [window_info.x_limits_plot(1), window_info.x_limits_plot(2)];

    fig = create_publishable_figure_local([120 90 1550 980]);
    tile = tiledlayout(fig, 3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('%s: subharmonic amplitude budget (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
        strrep(case_tag, '_', '\_'), sigma_value, t_index, t_value), ...
        'Interpreter', 'tex', 'FontSize', 15, 'FontWeight', 'bold');

    ax1 = nexttile(tile, 1);
    hold(ax1, 'on');
    plot(ax1, x_plot, ow3d.u_bare, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D bare u_{20}');
    plot(ax1, x_plot, mf12.u20, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, 'DisplayName', 'MF12 u_{20}');
    plot(ax1, x_plot, mf12.phix_bulk, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.5, 'DisplayName', 'MF12 (\phi_{20}^{bulk})_x');
    plot(ax1, x_plot, mf12.phix_surface, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 1.9, 'DisplayName', 'MF12 (\phi_{20}^{surf})_x');
    hold(ax1, 'off');
    style_axes_local(ax1, x_limits, compute_multi_series_ylimits_local({ow3d.u_bare, mf12.u20, mf12.phix_bulk, mf12.phix_surface}));
    xlabel(ax1, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax1, '$u_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax1, 'u_{20} overlay', 'Interpreter', 'tex');
    legend(ax1, 'Location', 'best', 'FontSize', 9);

    ax2 = nexttile(tile, 2);
    hold(ax2, 'on');
    plot(ax2, x_plot, ow3d.u_bare - mf12.u20, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, 'DisplayName', 'OW3D - MF12 u_{20}');
    plot(ax2, x_plot, ow3d.u_bare - mf12.phix_bulk, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.5, 'DisplayName', 'OW3D - MF12 (\phi_{20}^{bulk})_x');
    plot(ax2, x_plot, ow3d.u_bare - mf12.phix_surface, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 1.9, 'DisplayName', 'OW3D - MF12 (\phi_{20}^{surf})_x');
    hold(ax2, 'off');
    style_axes_local(ax2, x_limits, compute_multi_series_ylimits_local({ow3d.u_bare - mf12.u20, ow3d.u_bare - mf12.phix_bulk, ow3d.u_bare - mf12.phix_surface}));
    xlabel(ax2, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax2, 'Residual (m/s)', 'Interpreter', 'latex');
    title(ax2, 'u_{20} residuals', 'Interpreter', 'tex');
    legend(ax2, 'Location', 'best', 'FontSize', 9);

    ax3 = nexttile(tile, 3);
    hold(ax3, 'on');
    plot(ax3, x_plot, ow3d.u_bare, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D bare u_{20}');
    plot(ax3, x_plot, metrics.u20.mf12_u20.gain * mf12.u20, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, 'DisplayName', sprintf('g u_{20}, g=%.3f', metrics.u20.mf12_u20.gain));
    plot(ax3, x_plot, metrics.u20.bulk_phi_x.gain * mf12.phix_bulk, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.5, 'DisplayName', sprintf('g (\\phi_{20}^{bulk})_x, g=%.3f', metrics.u20.bulk_phi_x.gain));
    plot(ax3, x_plot, metrics.u20.surface_phi_x.gain * mf12.phix_surface, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 1.9, 'DisplayName', sprintf('g (\\phi_{20}^{surf})_x, g=%.3f', metrics.u20.surface_phi_x.gain));
    hold(ax3, 'off');
    style_axes_local(ax3, x_limits, compute_multi_series_ylimits_local({ow3d.u_bare, metrics.u20.mf12_u20.gain * mf12.u20, metrics.u20.bulk_phi_x.gain * mf12.phix_bulk, metrics.u20.surface_phi_x.gain * mf12.phix_surface}));
    xlabel(ax3, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax3, '$u_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax3, 'u_{20} gain-adjusted overlay', 'Interpreter', 'tex');
    legend(ax3, 'Location', 'best', 'FontSize', 9);

    ax4 = nexttile(tile, 4);
    hold(ax4, 'on');
    plot(ax4, x_plot, ow3d.w20, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D w_{20}');
    plot(ax4, x_plot, mf12.w20_z0, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, 'DisplayName', 'MF12 w_{20}(z=0)');
    plot(ax4, x_plot, mf12.w20_zeta11, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', 'MF12 w_{20}(z=\eta^{(1)})');
    hold(ax4, 'off');
    style_axes_local(ax4, x_limits, compute_multi_series_ylimits_local({ow3d.w20, mf12.w20_z0, mf12.w20_zeta11}));
    xlabel(ax4, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax4, '$w_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax4, 'w_{20} overlay', 'Interpreter', 'tex');
    legend(ax4, 'Location', 'best', 'FontSize', 9);

    ax5 = nexttile(tile, 5);
    hold(ax5, 'on');
    plot(ax5, x_plot, ow3d.w20 - mf12.w20_z0, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, 'DisplayName', 'OW3D - MF12 w_{20}(z=0)');
    plot(ax5, x_plot, ow3d.w20 - mf12.w20_zeta11, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', 'OW3D - MF12 w_{20}(z=\eta^{(1)})');
    hold(ax5, 'off');
    style_axes_local(ax5, x_limits, compute_multi_series_ylimits_local({ow3d.w20 - mf12.w20_z0, ow3d.w20 - mf12.w20_zeta11}));
    xlabel(ax5, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax5, 'Residual (m/s)', 'Interpreter', 'latex');
    title(ax5, 'w_{20} residuals', 'Interpreter', 'tex');
    legend(ax5, 'Location', 'best', 'FontSize', 9);

    ax6 = nexttile(tile, 6);
    hold(ax6, 'on');
    plot(ax6, x_plot, ow3d.w20, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D w_{20}');
    plot(ax6, x_plot, metrics.w20.mf12_w20_z0.gain * mf12.w20_z0, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.7, ...
        'DisplayName', sprintf('g w_{20}(z=0), g=%.3f', metrics.w20.mf12_w20_z0.gain));
    plot(ax6, x_plot, metrics.w20.mf12_w20_zeta11.gain * mf12.w20_zeta11, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, ...
        'DisplayName', sprintf('g w_{20}(z=\\eta^{(1)}), g=%.3f', metrics.w20.mf12_w20_zeta11.gain));
    hold(ax6, 'off');
    style_axes_local(ax6, x_limits, compute_multi_series_ylimits_local({ow3d.w20, metrics.w20.mf12_w20_z0.gain * mf12.w20_z0, metrics.w20.mf12_w20_zeta11.gain * mf12.w20_zeta11}));
    xlabel(ax6, '$x / \lambda$', 'Interpreter', 'latex'); ylabel(ax6, '$w_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax6, 'w_{20} gain-adjusted overlay', 'Interpreter', 'tex');
    legend(ax6, 'Location', 'best', 'FontSize', 9);

    annotation(fig, 'textbox', [0.11 0.005 0.87 0.06], ...
        'String', sprintf(['Local window = %.2f\\lambda around the |\\eta^{(1)}| envelope maximum. ' ...
        'Saved metrics include raw/gain-adjusted fits, local peak ratios, and w_{20} main-band peak ratios.'], ...
        window_info.window_lambda), ...
        'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontName', 'Times New Roman', 'FontSize', 10.5);

    exportgraphics(fig, fullfile(output_dir, ...
        sprintf('compare_ow3d_mf12_subharmonic_amplitude_budget_sigma_%03d_tidx_%04d.png', window_info.sigma_idx, t_index)), ...
        'Resolution', 300);
end

function style_axes_local(ax, x_limits, y_limits)
    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 11.5, 'LineWidth', 1.0);
    xlim(ax, x_limits);
    ylim(ax, y_limits);
end

function out = compute_mf12_subharmonic_targets_local(eta11, x_vec, depth, gravity, energy_keep)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    nx = numel(x_vec);
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

    kx = kx_grid(positive_idx).';
    ky = zeros(size(kx));
    a = 2 * real(eta_hat(positive_idx)).';
    b = 2 * imag(eta_hat(positive_idx)).';

    coeffs2 = mf12_direct_coefficients(2, gravity, depth, a, b, kx, ky, 0, 0, 0);
    [eta20, phi_surface, u20, ~, w20] = mf12_second_subharmonic_kinematics(coeffs2, x_vec.', 0, 0, 0);
    phi_bulk = mf12_second_subharmonic_bulk_phi_z0_local(coeffs2, x_vec.', 0, 0, 0);
    [u20_zeta11, w20_zeta11] = evaluate_mf12_subharmonic_pointwise_local(coeffs2, x_vec, eta11);

    out = struct();
    out.eta20 = eta20(:);
    out.coeffs2 = coeffs2;
    out.linear_indices = positive_idx(:);
    out.u20 = u20(:);
    out.phi_surface = phi_surface(:);
    out.phi_bulk = phi_bulk(:);
    out.phix_surface = spectral_derivative_x_local(out.phi_surface, x_vec);
    out.phix_bulk = spectral_derivative_x_local(out.phi_bulk, x_vec);
    out.w20 = w20(:);
    out.w20_z0 = w20(:);
    out.w20_zeta11 = w20_zeta11(:);
    out.u20_zeta11 = u20_zeta11(:);
end

function [u_zeta, w_zeta] = evaluate_mf12_subharmonic_pointwise_local(coeffs2, x_vec, z_vec)
    npts = numel(x_vec);
    u_zeta = zeros(npts, 1);
    w_zeta = zeros(npts, 1);
    for idx = 1:npts
        [~, ~, u_zeta(idx), ~, w_zeta(idx)] = mf12_second_subharmonic_kinematics( ...
            coeffs2, x_vec(idx), 0, z_vec(idx), 0);
    end
end

function phi = mf12_second_subharmonic_bulk_phi_z0_local(coeffs, x, y, z, t)
    phi = zeros(size(t));
    Z = z + coeffs.h;
    cnm = 0;
    for n = 1:coeffs.N
        for m = n+1:coeffs.N
            cnm = cnm + 1;
            idxMinus = 2 * cnm;
            theta_nm = coeffs.omega_npm(idxMinus) .* t ...
                - coeffs.kx_npm(idxMinus) .* x - coeffs.ky_npm(idxMinus) .* y;
            factorZ = coeffs.F_npm(idxMinus) .* cosh(coeffs.kappa_npm(idxMinus) .* Z);
            phi = phi + factorZ .* (coeffs.A_npm(idxMinus) .* sin(theta_nm) ...
                - coeffs.B_npm(idxMinus) .* cos(theta_nm));
        end
    end
end

function out = compute_ow3d_subharmonic_targets_local(eta_phases, vars_phases, h_vec, sigma_idx, sigma_value, x_vec, four_phase_coef, diff_alpha)
    eta_surface_phases = eta_phases;
    phi_surface_phases = squeeze(vars_phases.phi(:, sigma_idx, :));
    u_surface_phases = squeeze(vars_phases.u(:, sigma_idx, :));
    w_surface_phases = squeeze(vars_phases.w(:, sigma_idx, :));

    h_row = h_vec(:).';
    phix_sigma_phases = diffxeven_phasewise_local(phi_surface_phases, x_vec, diff_alpha);
    etax_phases = diffxeven_phasewise_local(eta_surface_phases, x_vec, diff_alpha);
    hx_row = diffxeven_local(h_row(:), x_vec, diff_alpha).';
    d_surface_phases = h_row + eta_surface_phases;

    chain_metric_phases = (((1 - sigma_value) .* hx_row ./ d_surface_phases) ...
        - sigma_value .* etax_phases ./ d_surface_phases) .* w_surface_phases;
    chain_surface_phases = -sigma_value .* etax_phases .* w_surface_phases;

    raw_harmonics = reconstruct_harmonics_1d_local(u_surface_phases, four_phase_coef);
    phix_harmonics = reconstruct_harmonics_1d_local(phix_sigma_phases, four_phase_coef);
    chain_metric_harmonics = reconstruct_harmonics_1d_local(chain_metric_phases, four_phase_coef);
    chain_surface_harmonics = reconstruct_harmonics_1d_local(chain_surface_phases, four_phase_coef);
    w_harmonics = reconstruct_harmonics_1d_local(w_surface_phases, four_phase_coef);

    closure_metric = raw_harmonics(4, :).'- (phix_harmonics(4, :).'+ chain_metric_harmonics(4, :).');
    closure_surface = raw_harmonics(4, :).'- (phix_harmonics(4, :).'+ chain_surface_harmonics(4, :).');
    if max(abs(closure_surface)) <= max(abs(closure_metric))
        selected_chain = chain_surface_harmonics(4, :).';
    else
        selected_chain = chain_metric_harmonics(4, :).';
    end

    out = struct();
    out.u_raw = raw_harmonics(4, :).';
    out.phix_sigma = phix_harmonics(4, :).';
    out.u_chain = selected_chain;
    out.u_bare = out.u_raw - out.u_chain;
    out.w20 = w_harmonics(4, :).';
end

function metric = compare_with_gain_local(reference, candidate, window_mask)
    reference = reference(:);
    candidate = candidate(:);
    metric.raw = compare_series_metrics_local(reference, candidate);
    denom = dot(candidate, candidate);
    if denom <= eps
        gain = 0;
    else
        gain = dot(reference, candidate) / denom;
    end
    candidate_adj = gain * candidate;
    metric.gain = gain;
    metric.gain_adjusted = compare_series_metrics_local(reference, candidate_adj);
    metric.local_peak_ratio = local_peak_ratio_local(reference, candidate, window_mask);
    metric.local_peak_ratio_gain_adjusted = local_peak_ratio_local(reference, candidate_adj, window_mask);
end

function ratio = local_peak_ratio_local(reference, candidate, mask)
    ref_local = reference(mask);
    cand_local = candidate(mask);
    ratio = max(abs(cand_local)) / max(abs(ref_local));
end

function ratio = compute_main_band_peak_ratio_local(reference, candidate, x_vec, kp, harmonic_order)
    [k_plot, amp_ref, amp_candidate] = compute_one_sided_spectrum_local(reference, candidate, x_vec, kp);
    band = harmonic_band_edges_local(harmonic_order, kp);
    mask = k_plot >= band(1) & k_plot <= band(2);
    if ~any(mask)
        ratio = NaN;
        return;
    end
    ratio = max(amp_candidate(mask)) / max(amp_ref(mask));
end

function [k_plot, amp_a, amp_b] = compute_one_sided_spectrum_local(field_a, field_b, x_vec, kp)
    field_a = field_a(:);
    field_b = field_b(:);
    x_vec = x_vec(:);
    nx = numel(field_a);
    dx = x_vec(2) - x_vec(1);
    fft_a = fft(field_a) / nx;
    fft_b = fft(field_b) / nx;
    if mod(nx, 2) == 0
        pos_idx = 1:(nx / 2 + 1);
    else
        pos_idx = 1:((nx + 1) / 2);
    end
    k_vals = vwa_kxgrid_local(nx, dx);
    k_plot = k_vals(pos_idx) / kp;
    amp_a = 2 * abs(fft_a(pos_idx));
    amp_b = 2 * abs(fft_b(pos_idx));
    amp_a(1) = abs(fft_a(1));
    amp_b(1) = abs(fft_b(1));
    if mod(nx, 2) == 0
        amp_a(end) = abs(fft_a(pos_idx(end)));
        amp_b(end) = abs(fft_b(pos_idx(end)));
    end
    k_plot = k_plot(:) * kp;
    amp_a = amp_a(:);
    amp_b = amp_b(:);
end

function info = select_local_window_local(eta11, x_vec, lambda, window_lambda, sigma_idx)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    envelope = abs(hilbert(eta11));
    [~, center_idx] = max(envelope);
    x_center = x_vec(center_idx);
    x_limits = [x_center - 0.5 * window_lambda * lambda, x_center + 0.5 * window_lambda * lambda];
    mask = x_vec >= x_limits(1) & x_vec <= x_limits(2);
    x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / lambda;
    info = struct();
    info.mask = mask;
    info.center_idx = center_idx;
    info.center_x = x_center;
    info.window_lambda = window_lambda;
    info.x_limits = x_limits;
    info.x_limits_plot = [x_plot(find(mask, 1, 'first')), x_plot(find(mask, 1, 'last'))];
    info.sigma_idx = sigma_idx;
end

function print_metric_line_local(label, metric)
    fprintf('%s : corr = %.6f, RMSE = %.6e, peak ratio = %.6f, gain = %.6f, gain-adjusted RMSE = %.6e, local peak = %.6f\n', ...
        label, metric.raw.corr, metric.raw.rmse, metric.raw.peak_ratio, metric.gain, ...
        metric.gain_adjusted.rmse, metric.local_peak_ratio);
end

function metrics = compare_series_metrics_local(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    metrics.rmse = sqrt(mean((reference - candidate).^2));
    corr_matrix = corrcoef(reference, candidate);
    metrics.corr = corr_matrix(1, 2);
    metrics.peak_ratio = max(abs(candidate)) / max(abs(reference));
end

function harmonics_out = filter_harmonics_x_only_local(harmonics_in, x_vec, kp)
    harmonics_out = harmonics_in;
    for n = 1:size(harmonics_in, 1)
        harmonics_out(n, :) = frequency_filtering_1d_local(squeeze(harmonics_in(n, :)), x_vec, kp, n);
    end
end

function field_out = frequency_filtering_1d_local(field_in, x_vec, kp, n)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    band_edges = harmonic_band_edges_local(n, kp);
    spectrum = fft(field_in);
    mask = (abs(kx) >= band_edges(1)) & (abs(kx) <= band_edges(2));
    field_out = real(ifft(spectrum .* mask));
end

function edges = harmonic_band_edges_local(order, kp)
    switch order
        case 1
            edges = [0.5, 1.5] * kp;
        case 2
            edges = [1.5, 2.5] * kp;
        case 3
            edges = [2.5, 3.5] * kp;
        case 4
            edges = [0.0, 0.75] * kp;
        otherwise
            error('Unsupported harmonic order %d.', order);
    end
end

function derivative = spectral_derivative_x_local(field_in, x_vec)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    kx = vwa_kxgrid_local(nx, dx);
    derivative = real(ifft(1i * kx .* fft(field_in)));
end

function derivative = diffxeven_phasewise_local(fields_by_phase, x_vec, alpha)
    derivative = zeros(size(fields_by_phase));
    for phase_idx = 1:size(fields_by_phase, 1)
        derivative(phase_idx, :) = diffxeven_local(fields_by_phase(phase_idx, :).', x_vec, alpha).';
    end
end

function derivative = diffxeven_local(field_in, x_vec, alpha)
    field_in = field_in(:);
    x_vec = x_vec(:);
    if nargin < 3 || isempty(alpha)
        alpha = 2;
    end
    nx = numel(field_in);
    if nx < 2
        derivative = zeros(size(field_in));
        return;
    end
    dx = diff(x_vec);
    alpha = max(1, min(alpha, floor((nx - 1) / 2)));
    rank = 2 * alpha + 1;
    coeff = build_stencil_even_local(alpha, 1) / dx(1);
    derivative = zeros(nx, 1);
    for ix = 1:alpha
        derivative(ix) = coeff(:, ix).' * field_in(1:rank);
    end
    for ix = alpha + 1:nx - alpha
        derivative(ix) = coeff(:, alpha + 1).' * field_in(ix - alpha:ix + alpha);
    end
    for ix = nx - alpha + 1:nx
        derivative(ix) = coeff(:, rank - (nx - ix)).' * field_in(nx - rank + 1:nx);
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

function harmonics = reconstruct_harmonics_1d_local(fields_by_phase, coef)
    analytic_part = hilbert(fields_by_phase.').';
    all_fields = cat(1, real(fields_by_phase), -imag(analytic_part));
    harmonics = zeros(4, size(fields_by_phase, 2));
    for n = 1:4
        harmonics(n, :) = coef(n, :) * all_fields;
    end
end

function data = read_ow3d_kinematics_snapshot_local(kin_path, phit_mode)
    [it, eta, ~, ~, phi, ~, ~, ~, u, ~, w, ~, ~, ~, x, y, h, sigma, t] = ...
        read_kinematics_file_local(kin_path, phit_mode); %#ok<ASGLU>
    data = struct('it', it, 'eta', eta, 'phi', phi, 'u', u, 'w', w, 'x', x, 'y', y, 'h', h, 'sigma', sigma, 't', t);
end

function [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = read_kinematics_file_local(file_path, phit_mode)
    nbits = 32;
    compute_derivatives = false;

    if nbits == 32
        int_nbit = 'int';
    else
        int_nbit = 'int64';
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
    h(:) = tmp(3:5:5 * nx * ny);

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

    eta = eta(1:it, :, :);
    phi = phi(1:it, :, :, :);
    u = u(1:it, :, :, :);
    v = v(1:it, :, :, :);
    w = w(1:it, :, :, :);
    uz = uz(1:it, :, :, :);
    vz = vz(1:it, :, :, :);
    wz = wz(1:it, :, :, :);
    t = t(1:it);

    etat_m = 0;
    etatt_m = 0;
    phit_m = 0;
    p_m = 0;
    ut_m = 0;
    if compute_derivatives
        error('Derivative reconstruction not used in this comparison script.');
    end

    switch lower(phit_mode)
        case 'uncorrected'
        case 'sigma_corrected'
        otherwise
            error('Unsupported phit_mode: %s', phit_mode);
    end
end

function idx = resolve_sigma_index_local(CFG, sigma_vec)
    if strcmpi(CFG.sigma_mode, 'index')
        idx = CFG.sigma_index;
    elseif strcmpi(CFG.sigma_mode, 'value')
        [~, idx] = min(abs(sigma_vec - CFG.sigma_value));
    else
        idx = numel(sigma_vec);
    end
end

function idx = resolve_time_index_local(CFG, requested_index, it_vec, nt)
    if ~isempty(requested_index)
        idx = requested_index;
    elseif ~isempty(CFG.time_index)
        idx = CFG.time_index;
    else
        idx = max(1, nt - CFG.default_time_index_from_end + 1);
    end

    if any(it_vec == idx)
        idx = find(it_vec == idx, 1, 'first');
    end

    idx = max(1, min(nt, idx));
end

function assert_phase_compatibility_local(phase_data, data_by_phase, idx)
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
    if ~isequal(size(phase_data.eta), size(ref.eta))
        error('Kinematics eta array size mismatch between phase 0 and phase index %d.', idx - 1);
    end
    if ~isequal(size(phase_data.phi), size(ref.phi)) || ~isequal(size(phase_data.u), size(ref.u))
        error('Kinematics field size mismatch between phase 0 and phase index %d.', idx - 1);
    end
end

function kin_path = resolve_kinematics_path_local(case_folder, file_id)
    kin_path = fullfile(case_folder, sprintf('Kinematics%02d.bin', file_id));
    if ~isfile(kin_path)
        error('Missing kinematics file: %s', kin_path);
    end
end

function kd = extract_kd_from_case_pattern_local(folder_pattern)
    token = regexp(folder_pattern, 'kd([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once');
    if isempty(token)
        error('Could not parse kd from folder pattern: %s', folder_pattern);
    end
    kd = str2double(token{1});
end

function case_tag = build_case_tag_local(folder_pattern)
    phi_token = regexp(folder_pattern, '^(.*)_phi_%d$', 'tokens', 'once');
    if isempty(phi_token)
        case_tag = regexprep(folder_pattern, '[\\/:*?""<>|]', '_');
    else
        case_tag = phi_token{1};
    end
end

function alpha = resolve_ow3d_diff_alpha_local(case_folder, default_alpha)
    if nargin < 2 || isempty(default_alpha)
        default_alpha = 2;
    end
    input_path = fullfile(case_folder, 'OceanWave3D.inp');
    alpha = default_alpha;
    if ~isfile(input_path)
        return;
    end
    text = fileread(input_path);
    token = regexp(text, 'alpha\s*=\s*([0-9]+)', 'tokens', 'once');
    if ~isempty(token)
        alpha = str2double(token{1});
    end
end

function limits = compute_multi_series_ylimits_local(series_cell)
    ymin = inf;
    ymax = -inf;
    for idx = 1:numel(series_cell)
        values = series_cell{idx};
        ymin = min(ymin, min(values(:)));
        ymax = max(ymax, max(values(:)));
    end
    padding = 0.08 * max(ymax - ymin, eps);
    limits = [ymin - padding, ymax + padding];
end

function fig = create_publishable_figure_local(position)
    fig = figure('Color', 'w', 'Position', position, 'Renderer', 'painters');
end

function kx = vwa_kxgrid_local(nx, dx)
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
end

function idx_out = select_energy_dominant_indices_local(eta_hat, candidate_idx, energy_keep)
    amplitudes = abs(eta_hat(candidate_idx)).^2;
    [~, sort_idx] = sort(amplitudes, 'descend');
    sorted_candidate_idx = candidate_idx(sort_idx);
    cumulative_energy = cumsum(amplitudes(sort_idx));
    total_energy = cumulative_energy(end);
    if total_energy <= 0
        idx_out = sorted_candidate_idx;
        return;
    end
    n_keep = find(cumulative_energy >= energy_keep * total_energy, 1, 'first');
    idx_out = sorted_candidate_idx(1:n_keep);
    idx_out = sort(idx_out(:), 'ascend');
end
