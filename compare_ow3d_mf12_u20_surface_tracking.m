% compare_ow3d_mf12_u20_surface_tracking.m
% Compare OW3D bare u20 against MF12 u20 evaluated at z = 0 and z = eta^(1)(x)
% over a local window centered on the linear envelope maximum.

clc;
clear;
close all;

CFG = struct();

CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check3');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
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
CFG.window_lambda = 6.0;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');

folder_pattern_env = getenv('OW3D_FOLDER_PATTERN');
if ~isempty(folder_pattern_env)
    CFG.folder_pattern = folder_pattern_env;
end

window_lambda_env = getenv('MF12_WINDOW_LAMBDA');
if ~isempty(window_lambda_env)
    value = str2double(window_lambda_env);
    if isfinite(value) && value > 0
        CFG.window_lambda = value;
    end
end

helper_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
if ~isfolder(helper_dir)
    error('Missing MF12 helper directory: %s', helper_dir);
end
addpath(helper_dir);

data_by_phase = cell(1, numel(CFG.phi_shifts_deg));
time_index_by_phase = zeros(1, numel(CFG.phi_shifts_deg));

for idx = 1:numel(CFG.phi_shifts_deg)
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
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
case_kd = extract_kd_from_case_pattern_local(CFG.folder_pattern);
depth_value = case_kd / CFG.kp_depth;
case_folder_ref = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(1)));
ow3d_diff_alpha = resolve_ow3d_diff_alpha_local(case_folder_ref, []);
case_tag = build_case_tag_local(CFG.folder_pattern);
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
kp = 2 * pi / CFG.lambda;
if CFG.apply_x_filter
    eta_harmonics = filter_harmonics_x_only_local(eta_harmonics, x_vec, kp);
end

sigma_idx = resolve_sigma_index_local(CFG, sigma_vec);
sigma_value = sigma_vec(sigma_idx);
eta11_surface = squeeze(eta_harmonics(1, :)).';
ow3d = compute_ow3d_second_subharmonic_surface_local( ...
    eta_phases, vars_phases, ref.h(:, 1), sigma_idx, sigma_value, x_vec, four_phase_coef, ow3d_diff_alpha);
mf12 = compute_mf12_second_subharmonic_surface_tracking_local( ...
    eta11_surface, x_vec, depth_value, CFG.gravity, CFG.linear_energy_keep);

[window_mask, x_window_limits, center_idx] = select_envelope_window_local(eta11_surface, x_vec, CFG.lambda, CFG.window_lambda);
x_win = x_vec(window_mask);
z0_win = zeros(size(x_win));
zeta11_win = eta11_surface(window_mask);
[~, ~, u0_vec_win] = mf12_second_subharmonic_kinematics(mf12.coeffs2, x_win.', 0, z0_win.', 0);
[~, ~, u_eta11_vec_win] = mf12_second_subharmonic_kinematics(mf12.coeffs2, x_win.', 0, zeta11_win.', 0);
[u0_win, u_eta11_win] = evaluate_mf12_u20_pointwise_local(mf12.coeffs2, x_win, z0_win, zeta11_win);

ow3d_bare_win = ow3d.u_bare(window_mask);
ow3d_raw_win = ow3d.u_raw(window_mask);
ow3d_chain_win = ow3d.u_chain(window_mask);
mf12_z0_win = u0_win(:);
mf12_zeta11_win = u_eta11_win(:);

metrics_z0 = compare_series_metrics_local(ow3d_bare_win, mf12_z0_win);
metrics_zeta11 = compare_series_metrics_local(ow3d_bare_win, mf12_zeta11_win);

fprintf('\n=== Local window u20 tracking diagnostics ===\n');
fprintf('Window center x = %.6f m, width = %.2f lambda, retained points = %d\n', ...
    x_vec(center_idx), CFG.window_lambda, numel(x_win));
fprintf('OW3D bare vs MF12 z=0     : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_z0.corr, metrics_z0.rmse, metrics_z0.peak_ratio);
fprintf('OW3D bare vs MF12 z=eta11 : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_zeta11.corr, metrics_zeta11.rmse, metrics_zeta11.peak_ratio);
fprintf('Pointwise-vs-vectorized consistency: max|z=0| = %.6e, max|z=eta11| = %.6e\n', ...
    max(abs(u0_win - u0_vec_win(:))), max(abs(u_eta11_win - u_eta11_vec_win(:))));

x_plot = (x_win - 0.5 * (x_vec(1) + x_vec(end))) / CFG.lambda;
y_limits_main = compute_pairwise_ylimits_local({ow3d_bare_win}, {mf12_z0_win; mf12_zeta11_win});
y_limits_decomp = compute_multi_series_ylimits_local({ow3d_raw_win, ow3d_chain_win, ow3d_bare_win, mf12_zeta11_win});

fig = create_publishable_figure_local([140 110 1500 920]);
tile = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile, sprintf('Local u20 comparison near envelope maximum (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
    sigma_value, selected_time_index, t_selected), ...
    'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

ax1 = nexttile(tile);
hold(ax1, 'on');
plot(ax1, x_plot, ow3d_bare_win, 'k-', 'LineWidth', 1.8, 'DisplayName', 'OW3D bare u_{20}');
plot(ax1, x_plot, mf12_z0_win, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12 u_{20}(z=0)');
plot(ax1, x_plot, mf12_zeta11_win, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'DisplayName', 'MF12 u_{20}(z=\eta^{(1)})');
hold(ax1, 'off');
grid(ax1, 'on'); box(ax1, 'on');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax1, [x_plot(1), x_plot(end)]);
ylim(ax1, y_limits_main(1, :));
xlabel(ax1, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax1, '$u_{20}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax1, 'OW3D bare vs MF12 tracking evaluation', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax1, 'Location', 'best', 'FontSize', 10);
text(ax1, 0.02, 0.95, sprintf('z=0: corr = %.3f, RMSE = %.2e\nz=\\eta^{(1)}: corr = %.3f, RMSE = %.2e', ...
    metrics_z0.corr, metrics_z0.rmse, metrics_zeta11.corr, metrics_zeta11.rmse), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

ax2 = nexttile(tile);
hold(ax2, 'on');
plot(ax2, x_plot, ow3d_raw_win, 'Color', [0.10 0.10 0.10], 'LineWidth', 1.6, 'DisplayName', 'OW3D raw u_{20}');
plot(ax2, x_plot, ow3d_chain_win, 'Color', [0.82 0.24 0.14], 'LineWidth', 1.6, 'DisplayName', 'OW3D chain u_{20}');
plot(ax2, x_plot, ow3d_bare_win, 'Color', [0.18 0.55 0.34], 'LineWidth', 1.6, 'DisplayName', 'OW3D bare u_{20}');
plot(ax2, x_plot, mf12_zeta11_win, '--', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'DisplayName', 'MF12 u_{20}(z=\eta^{(1)})');
hold(ax2, 'off');
grid(ax2, 'on'); box(ax2, 'on');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax2, [x_plot(1), x_plot(end)]);
ylim(ax2, y_limits_decomp);
xlabel(ax2, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax2, '$u_{20}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax2, 'OW3D raw/bare/chain decomposition in the local window', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax2, 'Location', 'best', 'FontSize', 10);

ax3 = nexttile(tile);
hold(ax3, 'on');
plot(ax3, x_plot, eta11_surface(window_mask), 'k-', 'LineWidth', 1.6, 'DisplayName', '\eta^{(1)}');
plot(ax3, x_plot, zeta11_win, '--', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', 'z=\eta^{(1)}');
hold(ax3, 'off');
grid(ax3, 'on'); box(ax3, 'on');
set(ax3, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax3, [x_plot(1), x_plot(end)]);
xlabel(ax3, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax3, '$z$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax3, 'Tracking height used for MF12 evaluation', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax3, 'Location', 'best', 'FontSize', 10);

annotation(fig, 'textbox', [0.12 0.01 0.84 0.05], ...
    'String', sprintf(['Window centered on the maximum of the $\\eta^{(1)}$ envelope; width = %.2f $\\lambda$. ' ...
    'MF12 uses the same extracted linear spectrum in both cases; only the evaluation height changes from $z=0$ to $z=\\eta^{(1)}(x)$.'], ...
    CFG.window_lambda), ...
    'Interpreter', 'latex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(fig, fullfile(output_dir, ...
    sprintf('compare_ow3d_mf12_u20_surface_tracking_sigma_%03d_tidx_%04d.png', sigma_idx, selected_time_index)), ...
    'Resolution', 300);

results = struct();
results.case_tag = case_tag;
results.t_index = selected_time_index;
results.t_value = t_selected;
results.sigma_idx = sigma_idx;
results.sigma_value = sigma_value;
results.window_lambda = CFG.window_lambda;
results.window_x_limits = x_window_limits;
results.eta11_window = eta11_surface(window_mask);
results.ow3d_bare_u20 = ow3d_bare_win;
results.ow3d_raw_u20 = ow3d_raw_win;
results.ow3d_chain_u20 = ow3d_chain_win;
results.mf12_u20_z0 = mf12_z0_win;
results.mf12_u20_zeta11 = mf12_zeta11_win;
results.mf12_u20_z0_vectorized = u0_vec_win(:);
results.mf12_u20_zeta11_vectorized = u_eta11_vec_win(:);
results.metrics_z0 = metrics_z0;
results.metrics_zeta11 = metrics_zeta11;
results.linear_indices = mf12.linear_indices;
save(fullfile(output_dir, sprintf('compare_ow3d_mf12_u20_surface_tracking_sigma_%03d_tidx_%04d.mat', sigma_idx, selected_time_index)), ...
    'results', '-v7.3');

disp('Local u20 surface-tracking comparison complete.');

function data = read_ow3d_kinematics_snapshot_local(kin_path, phit_mode)
    [it, eta, ~, ~, phi, ~, ~, ~, u, ~, w, ~, ~, ~, x, y, h, sigma, t] = ...
        read_kinematics_file_local(kin_path, phit_mode); %#ok<ASGLU>
    data = struct('it', it, 'eta', eta, 'phi', phi, 'u', u, 'w', w, 'x', x, 'y', y, 'h', h, 'sigma', sigma, 't', t);
end

function [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, h, sigma, t] = read_kinematics_file_local(file_path, phit_mode)
    nbits = 32;
    compute_derivatives = true;

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
        c = build_stencil_even_local_local(alpha, 1);
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
    if ~isequal(size(phase_data.u), size(ref.u))
        error('Kinematics velocity array size mismatch between phase 0 and phase index %d.', idx - 1);
    end
    if numel(phase_data.sigma) ~= numel(ref.sigma) || any(abs(phase_data.sigma(:) - ref.sigma(:)) > 1e-12)
        error('Sigma grid mismatch between phase 0 and phase index %d.', idx - 1);
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
    sigma_k = 0.5 * kp;
    mask = exp(-((abs(kx) - n * kp).^2) / (2 * sigma_k^2));
    field_out = ifft(fft(field_in) .* mask);
    if isreal(field_in)
        field_out = real(field_out);
    end
end

function sigma_idx = resolve_sigma_index_local(CFG, sigma_vec)
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

function out = compute_ow3d_second_subharmonic_surface_local(eta_phases, vars_phases, h_vec, sigma_idx, sigma_value, x_vec, four_phase_coef, diff_alpha)
    eta_surface_phases = eta_phases;
    phi_surface_phases = squeeze(vars_phases.phi(:, sigma_idx, :));
    u_surface_phases = squeeze(vars_phases.u(:, sigma_idx, :));
    w_surface_phases = squeeze(vars_phases.w(:, sigma_idx, :));
    h_row = h_vec(:).';

    phix_sigma_phases = diffxeven_phasewise_local(phi_surface_phases, x_vec, diff_alpha);
    etax_phases = diffxeven_phasewise_local(eta_surface_phases, x_vec, diff_alpha);
    hx_row = diffxeven_local(h_row(:), x_vec, diff_alpha).';

    d_surface_phases = max(h_row + eta_surface_phases, eps);
    hx_phases = repmat(hx_row, size(eta_surface_phases, 1), 1);
    chain_metric_phases = (((1 - sigma_value) .* hx_phases ./ d_surface_phases) ...
        - (sigma_value .* etax_phases ./ d_surface_phases)) .* w_surface_phases;
    chain_surface_phases = -sigma_value .* etax_phases .* w_surface_phases;

    raw_harmonics = reconstruct_harmonics_1d_local(u_surface_phases, four_phase_coef);
    phix_harmonics = reconstruct_harmonics_1d_local(phix_sigma_phases, four_phase_coef);
    chain_metric_harmonics = reconstruct_harmonics_1d_local(chain_metric_phases, four_phase_coef);
    chain_surface_harmonics = reconstruct_harmonics_1d_local(chain_surface_phases, four_phase_coef);

    closure_metric = raw_harmonics(4, :).'- (phix_harmonics(4, :).'+ chain_metric_harmonics(4, :).');
    closure_surface = raw_harmonics(4, :).'- (phix_harmonics(4, :).'+ chain_surface_harmonics(4, :).');
    if max(abs(closure_metric)) <= max(abs(closure_surface))
        selected_chain_harmonics = chain_metric_harmonics;
    else
        selected_chain_harmonics = chain_surface_harmonics;
    end

    out = struct();
    out.u_raw = raw_harmonics(4, :).';
    out.phix_sigma = phix_harmonics(4, :).';
    out.u_chain = selected_chain_harmonics(4, :).';
    out.u_bare = out.u_raw - out.u_chain;
end

function out = compute_mf12_second_subharmonic_surface_tracking_local(eta11, x_vec, depth, gravity, energy_keep)
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

    if isempty(positive_idx)
        error('No positive-k components selected for MF12.');
    end

    kx = kx_grid(positive_idx).';
    ky = zeros(size(kx));
    a = 2 * real(eta_hat(positive_idx)).';
    b = 2 * imag(eta_hat(positive_idx)).';
    coeffs2 = mf12_direct_coefficients(2, gravity, depth, a, b, kx, ky, 0, 0, 0);

    out = struct();
    out.coeffs2 = coeffs2;
    out.linear_indices = positive_idx(:);
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
    coeff = build_stencil_even_local_local(alpha, 1) / dx(1);
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

function fx = build_stencil_even_local_local(alpha, der)
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

function positive_idx = select_energy_dominant_indices_local(spectrum, candidate_idx, energy_keep)
    energy_keep = min(max(energy_keep, 0), 1);
    if isempty(candidate_idx)
        positive_idx = candidate_idx;
        return;
    end
    spectral_energy = abs(spectrum(candidate_idx)).^2;
    total_energy = sum(spectral_energy);
    if total_energy <= 0 || energy_keep >= 1
        positive_idx = candidate_idx;
        return;
    end
    [sorted_energy, order] = sort(spectral_energy, 'descend');
    cumulative_energy = cumsum(sorted_energy) / total_energy;
    cutoff_idx = find(cumulative_energy >= energy_keep, 1, 'first');
    positive_idx = sort(candidate_idx(order(1:cutoff_idx)));
end

function [u_z0, u_zeta] = evaluate_mf12_u20_pointwise_local(coeffs2, x_win, z0_win, zeta11_win)
    npts = numel(x_win);
    u_z0 = zeros(npts, 1);
    u_zeta = zeros(npts, 1);

    for idx = 1:npts
        [~, ~, u_z0(idx)] = mf12_second_subharmonic_kinematics(coeffs2, x_win(idx), 0, z0_win(idx), 0);
        [~, ~, u_zeta(idx)] = mf12_second_subharmonic_kinematics(coeffs2, x_win(idx), 0, zeta11_win(idx), 0);
    end
end

function [mask, x_limits, center_idx] = select_envelope_window_local(eta11, x_vec, lambda, window_lambda)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    envelope = abs(hilbert(eta11));
    [~, center_idx] = max(envelope);
    half_width = 0.5 * window_lambda * lambda;
    x_center = x_vec(center_idx);
    x_limits = [x_center - half_width, x_center + half_width];
    mask = x_vec >= x_limits(1) & x_vec <= x_limits(2);
end

function metrics = compare_series_metrics_local(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    metrics = struct('corr', NaN, 'rmse', NaN, 'peak_ratio', NaN);
    cc = corrcoef(reference, candidate);
    metrics.corr = cc(1, 2);
    metrics.rmse = sqrt(mean((reference - candidate).^2));
    ref_peak = max(abs(reference));
    if ref_peak > eps
        metrics.peak_ratio = max(abs(candidate)) / ref_peak;
    end
end

function y_limits = compute_pairwise_ylimits_local(fields_a, fields_b)
    values = [fields_a{1}(:); fields_b{1}(:); fields_b{2}(:)];
    y_abs_max = max(abs(values));
    if y_abs_max == 0
        y_abs_max = 1;
    end
    padding = 0.08 * y_abs_max;
    y_limits = [-y_abs_max - padding, y_abs_max + padding];
end

function y_limits = compute_multi_series_ylimits_local(series_list)
    values = [];
    for idx = 1:numel(series_list)
        values = [values; series_list{idx}(:)]; %#ok<AGROW>
    end
    max_abs = max(abs(values));
    if max_abs < eps
        y_limits = [-1, 1];
    else
        y_limits = 1.1 * [-max_abs, max_abs];
    end
end

function fig = create_publishable_figure_local(fig_position)
    fig = figure('Color', 'w', 'Position', fig_position, 'Renderer', 'painters');
end

function kin_path = resolve_kinematics_path_local(case_folder, file_id)
    if file_id < 10
        kin_path = fullfile(case_folder, sprintf('Kinematics0%d.bin', file_id));
    else
        kin_path = fullfile(case_folder, sprintf('Kinematics%d.bin', file_id));
    end
end

function time_index = resolve_time_index_local(CFG, cfg_time_index, n_times_valid, ~)
    if isempty(cfg_time_index)
        time_index = n_times_valid - CFG.default_time_index_from_end + 1;
        time_index = max(1, min(n_times_valid, time_index));
    elseif cfg_time_index > 0
        time_index = min(cfg_time_index, n_times_valid);
    else
        time_index = max(1, n_times_valid + cfg_time_index + 1);
    end
end

function alpha = resolve_ow3d_diff_alpha_local(case_folder, alpha_override)
    if nargin >= 2 && ~isempty(alpha_override)
        alpha = alpha_override;
        return;
    end
    inp_path = fullfile(case_folder, 'OceanWave3D.inp');
    alpha = 3;
    if ~isfile(inp_path)
        return;
    end
    lines = regexp(fileread(inp_path), '\r\n|\n|\r', 'split');
    lines = cellfun(@strtrim, lines, 'UniformOutput', false);
    lines = lines(~cellfun(@isempty, lines));
    if numel(lines) < 4
        return;
    end
    tokens = regexp(lines{4}, '[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?', 'match');
    if isempty(tokens)
        return;
    end
    alpha_value = str2double(tokens{1});
    if isfinite(alpha_value) && alpha_value >= 1
        alpha = round(alpha_value);
    end
end

function kd = extract_kd_from_case_pattern_local(folder_pattern)
    token = regexp(folder_pattern, 'kd(?<kd>\d+(?:\.\d+)?)', 'names', 'once');
    if isempty(token) || ~isfield(token, 'kd')
        error('Unable to parse kd from CFG.folder_pattern: %s', folder_pattern);
    end
    kd = str2double(token.kd);
end

function case_tag = build_case_tag_local(folder_pattern)
    case_tag = strrep(folder_pattern, '_phi_%d', '');
    case_tag = regexprep(case_tag, '[^A-Za-z0-9._-]', '_');
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
