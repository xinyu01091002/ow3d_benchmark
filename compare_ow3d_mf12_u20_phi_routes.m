% compare_ow3d_mf12_u20_phi_routes.m
% Compare OW3D bare u20 against MF12 u20, d/dx of MF12 bulk phi20(z=0),
% and d/dx of MF12 surface phi20.

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
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');

folder_pattern_env = getenv('OW3D_FOLDER_PATTERN');
if ~isempty(folder_pattern_env)
    CFG.folder_pattern = folder_pattern_env;
end

mf12_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
if ~isfolder(mf12_dir)
    error('Missing MF12 source directory: %s', mf12_dir);
end
addpath(mf12_dir);

case_kd = extract_kd_from_case_pattern_local(CFG.folder_pattern);
depth_value = case_kd / CFG.kp_depth;
case_folder_ref = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(1)));
case_tag = build_case_tag_local(CFG.folder_pattern);
output_dir = fullfile(CFG.output_dir, case_tag);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

source_mat = fullfile(output_dir, 'compare_ow3d_mf12_u20_from_phi_chain_sigma_018_tidx_0066.mat');
if ~isfile(source_mat)
    error('Missing prerequisite file: %s. Run compare_ow3d_mf12_u20_from_phi_chain.m first.', source_mat);
end
loaded = load(source_mat, 'results');
source = loaded.results;

kin_path = resolve_kinematics_path_local(case_folder_ref, CFG.kinematics_file_id);
x_vec = read_x_grid_local(kin_path);
selected_time_index = source.t_index;
t_selected = source.t_value;
sigma_idx = source.sigma_idx;
sigma_value = source.sigma_value;
ow3d = source.ow3d;
mf12 = compute_mf12_phi_routes_from_coeffs_local(source.mf12.coeffs2, x_vec);

fprintf('Loaded prerequisite diagnostics: %s\n', source_mat);
fprintf('Using kinematics time index %d (t = %.6f s)\n', selected_time_index, t_selected);
fprintf('Using depth h = %.6f m from kd = %.4f and kp = %.4f 1/m\n', ...
    depth_value, case_kd, CFG.kp_depth);

metrics_u20 = compare_series_metrics_local(ow3d.u_bare, mf12.u20);
metrics_bulk = compare_series_metrics_local(ow3d.u_bare, mf12.phix_bulk);
metrics_surface = compare_series_metrics_local(ow3d.u_bare, mf12.phix_surface);
metrics_surface_vs_u20 = compare_series_metrics_local(mf12.u20, mf12.phix_surface);
metrics_bulk_vs_u20 = compare_series_metrics_local(mf12.u20, mf12.phix_bulk);

fprintf('\n=== MF12 phi-route diagnostics ===\n');
fprintf('OW3D bare vs MF12 u20            : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_u20.corr, metrics_u20.rmse, metrics_u20.peak_ratio);
fprintf('OW3D bare vs MF12 bulk phi20_x   : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_bulk.corr, metrics_bulk.rmse, metrics_bulk.peak_ratio);
fprintf('OW3D bare vs MF12 surface phi20_x: corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_surface.corr, metrics_surface.rmse, metrics_surface.peak_ratio);
fprintf('MF12 surface phi20_x vs MF12 u20 : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_surface_vs_u20.corr, metrics_surface_vs_u20.rmse, metrics_surface_vs_u20.peak_ratio);
fprintf('MF12 bulk phi20_x vs MF12 u20    : corr = %.6f, RMSE = %.6e, peak ratio = %.6f\n', ...
    metrics_bulk_vs_u20.corr, metrics_bulk_vs_u20.rmse, metrics_bulk_vs_u20.peak_ratio);

x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / CFG.lambda;
x_limits = resolve_plot_xlim_local(x_plot, ow3d.u_bare, 5.0);
y_limits_main = compute_multi_series_ylimits_local({ow3d.u_bare, mf12.u20, mf12.phix_bulk, mf12.phix_surface});
y_limits_diff = compute_multi_series_ylimits_local({mf12.phix_surface - mf12.u20, mf12.phix_bulk - mf12.u20});

fig = create_publishable_figure_local([140 110 1550 980]);
tile = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile, sprintf('OW3D bare u_{20} vs MF12 u_{20}/\\phi_{20,x} routes (\\sigma = %.3f, t index = %d, t = %.4f s)', ...
    sigma_value, selected_time_index, t_selected), ...
    'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

ax1 = nexttile(tile);
hold(ax1, 'on');
plot(ax1, x_plot, ow3d.u_bare, 'k-', 'LineWidth', 1.9, 'DisplayName', 'OW3D bare u_{20}');
plot(ax1, x_plot, mf12.u20, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'MF12 u_{20}');
plot(ax1, x_plot, mf12.phix_bulk, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', 'MF12 (\phi_{20}^{bulk})_x');
plot(ax1, x_plot, mf12.phix_surface, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 2.0, 'DisplayName', 'MF12 (\phi_{20}^{surf})_x');
hold(ax1, 'off');
grid(ax1, 'on'); box(ax1, 'on');
set(ax1, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax1, x_limits);
ylim(ax1, y_limits_main);
xlabel(ax1, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax1, '$u_{20}$ (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax1, 'Direct comparison against OW3D bare u_{20}', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax1, 'Location', 'best', 'FontSize', 10);

ax2 = nexttile(tile);
hold(ax2, 'on');
plot(ax2, x_plot, mf12.phix_bulk - mf12.u20, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8, 'DisplayName', '(\phi_{20}^{bulk})_x - u_{20}');
plot(ax2, x_plot, mf12.phix_surface - mf12.u20, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 2.0, 'DisplayName', '(\phi_{20}^{surf})_x - u_{20}');
hold(ax2, 'off');
grid(ax2, 'on'); box(ax2, 'on');
set(ax2, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax2, x_limits);
ylim(ax2, y_limits_diff);
xlabel(ax2, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax2, 'Difference (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax2, 'Internal MF12 route differences', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax2, 'Location', 'best', 'FontSize', 10);

ax3 = nexttile(tile);
hold(ax3, 'on');
plot(ax3, x_plot, ow3d.u_bare - mf12.u20, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8, 'DisplayName', 'OW3D bare - MF12 u_{20}');
plot(ax3, x_plot, ow3d.u_bare - mf12.phix_bulk, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', 'OW3D bare - MF12 (\phi_{20}^{bulk})_x');
plot(ax3, x_plot, ow3d.u_bare - mf12.phix_surface, ':', 'Color', [0.18 0.55 0.34], 'LineWidth', 2.0, 'DisplayName', 'OW3D bare - MF12 (\phi_{20}^{surf})_x');
hold(ax3, 'off');
grid(ax3, 'on'); box(ax3, 'on');
set(ax3, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
xlim(ax3, x_limits);
ylim(ax3, compute_multi_series_ylimits_local({ow3d.u_bare - mf12.u20, ow3d.u_bare - mf12.phix_bulk, ow3d.u_bare - mf12.phix_surface}));
xlabel(ax3, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 13);
ylabel(ax3, 'Difference (m/s)', 'Interpreter', 'latex', 'FontSize', 13);
title(ax3, 'Residuals against OW3D bare u_{20}', 'Interpreter', 'tex', 'FontSize', 13, 'FontWeight', 'normal');
legend(ax3, 'Location', 'best', 'FontSize', 10);

annotation(fig, 'textbox', [0.12 0.01 0.84 0.05], ...
    'String', ['MF12 bulk \phi_{20} is reconstructed from the F_{n-m}cosh(\kappa h) branch before the z=0 surface-potential overwrite. ' ...
    'This separates the route mismatch from the MF12-vs-OW3D mismatch.'], ...
    'Interpreter', 'tex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(fig, fullfile(output_dir, ...
    sprintf('compare_ow3d_mf12_u20_phi_routes_sigma_%03d_tidx_%04d.png', sigma_idx, selected_time_index)), ...
    'Resolution', 300);

results = struct();
results.case_tag = case_tag;
results.t_index = selected_time_index;
results.t_value = t_selected;
results.sigma_idx = sigma_idx;
results.sigma_value = sigma_value;
results.ow3d = ow3d;
results.mf12 = mf12;
results.metrics_u20 = metrics_u20;
results.metrics_bulk = metrics_bulk;
results.metrics_surface = metrics_surface;
results.metrics_surface_vs_u20 = metrics_surface_vs_u20;
results.metrics_bulk_vs_u20 = metrics_bulk_vs_u20;
save(fullfile(output_dir, sprintf('compare_ow3d_mf12_u20_phi_routes_sigma_%03d_tidx_%04d.mat', sigma_idx, selected_time_index)), ...
    'results', '-v7.3');

disp('MF12 phi-route u20 comparison complete.');

function out = compute_mf12_phi_routes_from_coeffs_local(coeffs2, x_vec)
    x_vec = x_vec(:);
    [eta20, phi_surface, u20, ~, w20] = mf12_second_subharmonic_kinematics(coeffs2, x_vec.', 0, 0, 0);
    phi_bulk = mf12_second_subharmonic_bulk_phi_z0_local(coeffs2, x_vec.', 0, 0, 0);

    out = struct();
    out.eta20 = eta20(:);
    out.u20 = u20(:);
    out.w20 = w20(:);
    out.phi_surface = phi_surface(:);
    out.phi_bulk = phi_bulk(:);
    out.phix_surface = spectral_derivative_x_local(out.phi_surface, x_vec);
    out.phix_bulk = spectral_derivative_x_local(out.phi_bulk, x_vec);
    out.coeffs2 = coeffs2;
end

function x = read_x_grid_local(file_path)
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
    fread(fid, 1, 'int');
    fread(fid, 1, 'int');
    fread(fid, 1, 'int');
    fread(fid, 1, 'double');
    fread(fid, 1, 'int');
    fread(fid, 2, 'int');

    nx = floor((xend - xbeg) / xstride) + 1;
    ny = floor((yend - ybeg) / ystride) + 1;
    tmp = fread(fid, 5 * nx * ny, 'double');
    x_grid = zeros(nx, ny);
    x_grid(:) = tmp(1:5:5 * nx * ny);
    x = x_grid(:, 1);
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

function derivative = spectral_derivative_x_local(field_in, x_vec)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    kx = vwa_kxgrid_local(nx, dx);
    derivative = real(ifft(1i * kx .* fft(field_in)));
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

function out = compute_ow3d_second_subharmonic_surface_local(eta_phases, vars_phases, h_vec, sigma_idx, sigma_value, x_vec, four_phase_coef, diff_alpha)
    eta_surface_phases = eta_phases;
    phi_surface_phases = squeeze(vars_phases.phi(:, sigma_idx, :));
    u_surface_phases = squeeze(vars_phases.u(:, sigma_idx, :));
    w_surface_phases = squeeze(vars_phases.w(:, sigma_idx, :));

    h_row = h_vec(:).';
    phix_sigma_phases = diffxeven_phasewise_local(phi_surface_phases, x_vec, diff_alpha);
    etax_phases = diffxeven_phasewise_local(eta_surface_phases, x_vec, diff_alpha);
    hx_row = diffxeven_local(h_row(:), x_vec, diff_alpha).';

    d_surface_phases = h_row + eta_surface_phases;
    metric_metric = ((1 - sigma_value) .* hx_row - sigma_value .* etax_phases) ./ d_surface_phases;
    metric_surface = -sigma_value .* etax_phases;
    chain_metric_phases = metric_metric .* w_surface_phases;
    chain_surface_phases = metric_surface .* w_surface_phases;

    raw_harmonics = reconstruct_harmonics_1d_local(u_surface_phases, four_phase_coef);
    phix_harmonics = reconstruct_harmonics_1d_local(phix_sigma_phases, four_phase_coef);
    chain_metric_harmonics = reconstruct_harmonics_1d_local(chain_metric_phases, four_phase_coef);
    chain_surface_harmonics = reconstruct_harmonics_1d_local(chain_surface_phases, four_phase_coef);

    residual_metric = max(abs(raw_harmonics(4, :) - (phix_harmonics(4, :) + chain_metric_harmonics(4, :))));
    residual_surface = max(abs(raw_harmonics(4, :) - (phix_harmonics(4, :) + chain_surface_harmonics(4, :))));

    if residual_surface <= residual_metric
        selected_chain = chain_surface_harmonics(4, :).';
    else
        selected_chain = chain_metric_harmonics(4, :).';
    end

    out = struct();
    out.u_raw = raw_harmonics(4, :).';
    out.phix_sigma = phix_harmonics(4, :).';
    out.u_chain = selected_chain;
    out.u_bare = out.u_raw - out.u_chain;
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
    derivative = zeros(size(field_in));
    if nx < 3
        return;
    end
    dx = diff(x_vec);
    alpha = max(1, min(alpha, floor((nx - 1) / 2)));
    rank = 2 * alpha + 1;

    interior_offsets = -alpha:alpha;
    interior_weights = build_stencil_even_local(interior_offsets);
    for idx = (alpha + 1):(nx - alpha)
        span = idx + interior_offsets;
        derivative(idx) = interior_weights * field_in(span) / mean(dx(span(1:end-1)));
    end

    for idx = 1:alpha
        left_offsets = (1:rank) - idx;
        left_weights = build_stencil_even_local(left_offsets);
        derivative(idx) = left_weights * field_in(idx + left_offsets) / mean(dx(idx:(idx + rank - 2)));

        right_offsets = -fliplr(1:rank) + (nx - idx + 1);
        right_weights = build_stencil_even_local(right_offsets);
        center = nx - idx + 1;
        derivative(center) = right_weights * field_in(center + right_offsets) / mean(dx((center - rank + 1):(center - 1)));
    end
end

function weights = build_stencil_even_local(offsets)
    offsets = offsets(:);
    n = numel(offsets);
    A = zeros(n, n);
    b = zeros(n, 1);
    for row = 1:n
        A(row, :) = (offsets.').^(row - 1);
    end
    b(2) = 1;
    weights = (A \ b).';
end

function harmonics = reconstruct_harmonics_1d_local(phase_fields, coef)
    phase_fields = phase_fields.';
    analytic_signal = hilbert(phase_fields);
    phase_hilbert = imag(analytic_signal).';
    base_matrix = [phase_fields, phase_hilbert];
    harmonics = coef * base_matrix;
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
        idx = max(1, nt - CFG.default_time_index_from_end);
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
    if ~isequal(size(phase_data.eta), size(ref.eta)) || ~isequal(size(phase_data.phi), size(ref.phi))
        error('OW3D phase files are not size-compatible.');
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

function metrics = compare_series_metrics_local(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    metrics.rmse = sqrt(mean((reference - candidate).^2));
    corr_matrix = corrcoef(reference, candidate);
    metrics.corr = corr_matrix(1, 2);
    metrics.peak_ratio = max(abs(candidate)) / max(abs(reference));
end

function x_limits = resolve_plot_xlim_local(x_plot, eta11, half_width_lambda)
    eta11 = eta11(:);
    envelope = abs(hilbert(eta11));
    [~, idx_max] = max(envelope);
    x_center = x_plot(idx_max);
    x_limits = [x_center - half_width_lambda, x_center + half_width_lambda];
    x_limits(1) = max(x_limits(1), x_plot(1));
    x_limits(2) = min(x_limits(2), x_plot(end));
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
