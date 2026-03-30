% compare_ow3d_mf12_w20_sigma_scan.m
% Scan OW3D w20 across sigma levels and compare against MF12 w20 evaluated
% at z = 0 and z = eta^(1)(x).

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
CFG.linear_energy_keep = 0.99999;
CFG.top_sigma_layers = 6;
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
    process_case_local(CFG, CFG.case_patterns{case_idx});
end

disp('W20 sigma scan complete.');

function process_case_local(CFG, folder_pattern)
    data_by_phase = cell(1, numel(CFG.phi_shifts_deg));
    time_index_by_phase = zeros(1, numel(CFG.phi_shifts_deg));

    for idx = 1:numel(CFG.phi_shifts_deg)
        case_folder = fullfile(CFG.data_root, sprintf(folder_pattern, CFG.phi_shifts_deg(idx)));
        kin_path = resolve_kinematics_path_local(case_folder, CFG.kinematics_file_id);
        phase_data = read_ow3d_kinematics_snapshot_local(kin_path, CFG.phit_mode);
        assert_phase_compatibility_local(phase_data, data_by_phase, idx);
        data_by_phase{idx} = phase_data;
        time_index_by_phase(idx) = resolve_time_index_local(CFG, phase_data.it, numel(phase_data.t));
    end

    if any(time_index_by_phase ~= time_index_by_phase(1))
        error('Resolved time indices are inconsistent across the four phase files.');
    end

    selected_time_index = time_index_by_phase(1);
    ref = data_by_phase{1};
    x_vec = ref.x(:, 1);
    sigma_vec = ref.sigma(:);
    t_selected = ref.t(selected_time_index);
    case_kd = extract_kd_from_case_pattern_local(folder_pattern);
    depth_value = case_kd / CFG.kp_depth;
    case_tag = build_case_tag_local(folder_pattern);
    output_dir = fullfile(CFG.output_dir, case_tag);
    if ~isfolder(output_dir)
        mkdir(output_dir);
    end

    four_phase_coef = [
        0.25  0    -0.25  0     0    -0.25  0     0.25;
        0.25 -0.25  0.25 -0.25  0     0      0     0;
        0.25  0    -0.25  0     0     0.25  0    -0.25;
        0.25  0.25  0.25  0.25  0     0      0     0];

    eta_phases = zeros(numel(CFG.phi_shifts_deg), size(ref.eta, 2));
    w_phases = zeros(numel(CFG.phi_shifts_deg), numel(sigma_vec), size(ref.w, 3));
    for idx = 1:numel(CFG.phi_shifts_deg)
        phase_data = data_by_phase{idx};
        eta_phases(idx, :) = squeeze(phase_data.eta(selected_time_index, :, 1));
        w_phases(idx, :, :) = squeeze(phase_data.w(selected_time_index, :, :, 1));
    end

    eta_harmonics = reconstruct_harmonics_1d_local(eta_phases, four_phase_coef);
    kp = 2 * pi / CFG.lambda;
    if CFG.apply_x_filter
        eta_harmonics = filter_harmonics_x_only_local(eta_harmonics, x_vec, kp);
    end
    eta11_surface = eta_harmonics(1, :).';

    mf12 = compute_mf12_w20_local(eta11_surface, x_vec, depth_value, CFG.gravity, CFG.linear_energy_keep);

    nz = numel(sigma_vec);
    sigma_start = max(1, nz - CFG.top_sigma_layers + 1);
    sigma_indices = sigma_start:nz;

    scan = struct([]);
    for i = 1:numel(sigma_indices)
        sigma_idx = sigma_indices(i);
        w_surface_phases = squeeze(w_phases(:, sigma_idx, :));
        w_harmonics = reconstruct_harmonics_1d_local(w_surface_phases, four_phase_coef);
        ow_w20 = w_harmonics(4, :).';
        scan(i).sigma_idx = sigma_idx;
        scan(i).sigma_value = sigma_vec(sigma_idx);
        scan(i).ow3d_w20 = ow_w20;
        z_mean_sigma = depth_value * (scan(i).sigma_value - 1);
        z_moving_sigma = z_mean_sigma + scan(i).sigma_value * eta11_surface;
        mf12_w_sigma_mean = evaluate_mf12_w20_at_z_local(mf12.coeffs2, x_vec, z_mean_sigma + 0 * x_vec);
        mf12_w_sigma_moving = evaluate_mf12_w20_at_z_local(mf12.coeffs2, x_vec, z_moving_sigma);

        scan(i).z_mean_sigma = z_mean_sigma;
        scan(i).metrics_z0 = compare_series_metrics_local(ow_w20, mf12.w20_z0);
        scan(i).metrics_zeta11 = compare_series_metrics_local(ow_w20, mf12.w20_zeta11);
        scan(i).metrics_sigma_mean = compare_series_metrics_local(ow_w20, mf12_w_sigma_mean);
        scan(i).metrics_sigma_moving = compare_series_metrics_local(ow_w20, mf12_w_sigma_moving);
        scan(i).gain_z0 = compute_gain_local(ow_w20, mf12.w20_z0);
        scan(i).gain_zeta11 = compute_gain_local(ow_w20, mf12.w20_zeta11);
        scan(i).gain_sigma_mean = compute_gain_local(ow_w20, mf12_w_sigma_mean);
        scan(i).gain_sigma_moving = compute_gain_local(ow_w20, mf12_w_sigma_moving);
    end

    fprintf('\n=== W20 sigma scan: %s ===\n', case_tag);
    for i = 1:numel(scan)
        fprintf(['sigma[%02d]=%.4f z_mean=%.3f | z=0 corr=%.6f peak=%.6f gain=%.6f' ...
            ' | z=eta11 corr=%.6f peak=%.6f gain=%.6f' ...
            ' | z=sigma*h corr=%.6f peak=%.6f gain=%.6f' ...
            ' | z=sigma(h+eta)-h corr=%.6f peak=%.6f gain=%.6f\n'], ...
            scan(i).sigma_idx, scan(i).sigma_value, scan(i).z_mean_sigma, ...
            scan(i).metrics_z0.corr, scan(i).metrics_z0.peak_ratio, scan(i).gain_z0, ...
            scan(i).metrics_zeta11.corr, scan(i).metrics_zeta11.peak_ratio, scan(i).gain_zeta11, ...
            scan(i).metrics_sigma_mean.corr, scan(i).metrics_sigma_mean.peak_ratio, scan(i).gain_sigma_mean, ...
            scan(i).metrics_sigma_moving.corr, scan(i).metrics_sigma_moving.peak_ratio, scan(i).gain_sigma_moving);
    end

    create_sigma_scan_figure_local(output_dir, case_tag, t_selected, selected_time_index, x_vec, eta11_surface, scan, mf12);

    results = struct();
    results.case_tag = case_tag;
    results.t_index = selected_time_index;
    results.t_value = t_selected;
    results.mf12 = mf12;
    results.scan = scan;
    save(fullfile(output_dir, sprintf('compare_ow3d_mf12_w20_sigma_scan_tidx_%04d.mat', selected_time_index)), ...
        'results', '-v7.3');
end

function create_sigma_scan_figure_local(output_dir, case_tag, t_value, t_index, x_vec, eta11_surface, scan, mf12)
    sigma_values = arrayfun(@(s) s.sigma_value, scan);
    corr_z0 = arrayfun(@(s) s.metrics_z0.corr, scan);
    corr_zeta = arrayfun(@(s) s.metrics_zeta11.corr, scan);
    corr_sigma_mean = arrayfun(@(s) s.metrics_sigma_mean.corr, scan);
    corr_sigma_moving = arrayfun(@(s) s.metrics_sigma_moving.corr, scan);
    peak_z0 = arrayfun(@(s) s.metrics_z0.peak_ratio, scan);
    peak_zeta = arrayfun(@(s) s.metrics_zeta11.peak_ratio, scan);
    peak_sigma_mean = arrayfun(@(s) s.metrics_sigma_mean.peak_ratio, scan);
    peak_sigma_moving = arrayfun(@(s) s.metrics_sigma_moving.peak_ratio, scan);
    gain_z0 = arrayfun(@(s) s.gain_z0, scan);
    gain_zeta = arrayfun(@(s) s.gain_zeta11, scan);
    gain_sigma_mean = arrayfun(@(s) s.gain_sigma_mean, scan);
    gain_sigma_moving = arrayfun(@(s) s.gain_sigma_moving, scan);

    best_idx = numel(scan);
    x_ref = x_vec(:);
    w_sigma_mean_best = evaluate_mf12_w20_at_z_local(mf12.coeffs2, x_ref, scan(best_idx).z_mean_sigma + 0 * x_ref);
    w_sigma_moving_best = evaluate_mf12_w20_at_z_local(mf12.coeffs2, x_ref, scan(best_idx).sigma_value * eta11_surface(:) + scan(best_idx).z_mean_sigma);

    fig = figure('Color', 'w', 'Position', [120 100 1400 900], 'Renderer', 'painters');
    tile = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('%s: OW3D/MF12 w_{20} sigma scan (t index = %d, t = %.4f s)', ...
        strrep(case_tag, '_', '\_'), t_index, t_value), 'Interpreter', 'tex', 'FontSize', 15, 'FontWeight', 'bold');

    ax1 = nexttile(tile, 1);
    hold(ax1, 'on');
    plot(ax1, sigma_values, corr_z0, 'o-', 'LineWidth', 1.6, 'DisplayName', 'corr vs MF12 w_{20}(z=0)');
    plot(ax1, sigma_values, corr_zeta, 's-', 'LineWidth', 1.6, 'DisplayName', 'corr vs MF12 w_{20}(z=\eta^{(1)})');
    plot(ax1, sigma_values, corr_sigma_mean, 'd-', 'LineWidth', 1.6, 'DisplayName', 'corr vs MF12 w_{20}(z=h(\sigma-1))');
    plot(ax1, sigma_values, corr_sigma_moving, '^-', 'LineWidth', 1.6, 'DisplayName', 'corr vs MF12 w_{20}(z=\sigma(h+\eta^{(1)})-h)');
    hold(ax1, 'off');
    grid(ax1, 'on'); box(ax1, 'on');
    xlabel(ax1, '\sigma'); ylabel(ax1, 'corr');
    title(ax1, 'Correlation by sigma layer');
    legend(ax1, 'Location', 'best');

    ax2 = nexttile(tile, 2);
    hold(ax2, 'on');
    plot(ax2, sigma_values, peak_z0, 'o-', 'LineWidth', 1.6, 'DisplayName', 'peak ratio vs z=0');
    plot(ax2, sigma_values, peak_zeta, 's-', 'LineWidth', 1.6, 'DisplayName', 'peak ratio vs z=\eta^{(1)}');
    plot(ax2, sigma_values, peak_sigma_mean, 'd-', 'LineWidth', 1.6, 'DisplayName', 'peak ratio vs z=h(\sigma-1)');
    plot(ax2, sigma_values, peak_sigma_moving, '^-', 'LineWidth', 1.6, 'DisplayName', 'peak ratio vs z=\sigma(h+\eta^{(1)})-h');
    plot(ax2, sigma_values, gain_z0, 'o--', 'LineWidth', 1.3, 'DisplayName', 'gain vs z=0');
    plot(ax2, sigma_values, gain_zeta, 's--', 'LineWidth', 1.3, 'DisplayName', 'gain vs z=\eta^{(1)}');
    plot(ax2, sigma_values, gain_sigma_mean, 'd--', 'LineWidth', 1.3, 'DisplayName', 'gain vs z=h(\sigma-1)');
    plot(ax2, sigma_values, gain_sigma_moving, '^--', 'LineWidth', 1.3, 'DisplayName', 'gain vs z=\sigma(h+\eta^{(1)})-h');
    hold(ax2, 'off');
    grid(ax2, 'on'); box(ax2, 'on');
    xlabel(ax2, '\sigma'); ylabel(ax2, 'amplitude ratio / gain');
    title(ax2, 'Amplitude by sigma layer');
    legend(ax2, 'Location', 'best');

    ax3 = nexttile(tile, 3);
    hold(ax3, 'on');
    plot(ax3, x_ref, scan(best_idx).ow3d_w20, 'k-', 'LineWidth', 1.8, 'DisplayName', sprintf('OW3D w_{20}, \\sigma=%.4f', scan(best_idx).sigma_value));
    plot(ax3, x_ref, mf12.w20_z0, '--', 'LineWidth', 1.6, 'DisplayName', 'MF12 w_{20}(z=0)');
    plot(ax3, x_ref, mf12.w20_zeta11, '-', 'LineWidth', 1.6, 'DisplayName', 'MF12 w_{20}(z=\eta^{(1)})');
    plot(ax3, x_ref, w_sigma_mean_best, ':', 'LineWidth', 1.6, 'DisplayName', 'MF12 w_{20}(z=h(\sigma-1))');
    plot(ax3, x_ref, w_sigma_moving_best, '-.', 'LineWidth', 1.6, 'DisplayName', 'MF12 w_{20}(z=\sigma(h+\eta^{(1)})-h)');
    hold(ax3, 'off');
    grid(ax3, 'on'); box(ax3, 'on');
    xlabel(ax3, 'x (m)'); ylabel(ax3, 'w_{20} (m/s)');
    title(ax3, 'Top sigma layer overlay');
    legend(ax3, 'Location', 'best');

    ax4 = nexttile(tile, 4);
    hold(ax4, 'on');
    for i = 1:numel(scan)
        plot(ax4, x_ref, scan(i).ow3d_w20, 'LineWidth', 1.1, 'DisplayName', sprintf('\\sigma=%.4f', scan(i).sigma_value));
    end
    hold(ax4, 'off');
    grid(ax4, 'on'); box(ax4, 'on');
    xlabel(ax4, 'x (m)'); ylabel(ax4, 'OW3D w_{20} (m/s)');
    title(ax4, 'OW3D w_{20} across sigma layers');
    legend(ax4, 'Location', 'best');

    exportgraphics(fig, fullfile(output_dir, sprintf('compare_ow3d_mf12_w20_sigma_scan_tidx_%04d.png', t_index)), 'Resolution', 300);
end

function out = compute_mf12_w20_local(eta11, x_vec, depth, gravity, energy_keep)
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
    [~, ~, ~, ~, w20_z0] = mf12_second_subharmonic_kinematics(coeffs2, x_vec.', 0, 0, 0);
    [~, w20_zeta11] = evaluate_mf12_w20_pointwise_local(coeffs2, x_vec, eta11);

    out = struct();
    out.coeffs2 = coeffs2;
    out.w20_z0 = w20_z0(:);
    out.w20_zeta11 = w20_zeta11(:);
end

function [w_z0, w_zeta] = evaluate_mf12_w20_pointwise_local(coeffs2, x_vec, z_vec)
    npts = numel(x_vec);
    w_z0 = zeros(npts, 1);
    w_zeta = zeros(npts, 1);
    for idx = 1:npts
        [~, ~, ~, ~, w_z0(idx)] = mf12_second_subharmonic_kinematics(coeffs2, x_vec(idx), 0, 0, 0);
        [~, ~, ~, ~, w_zeta(idx)] = mf12_second_subharmonic_kinematics(coeffs2, x_vec(idx), 0, z_vec(idx), 0);
    end
end

function w_out = evaluate_mf12_w20_at_z_local(coeffs2, x_vec, z_vec)
    [~, ~, ~, ~, w_out] = mf12_second_subharmonic_kinematics(coeffs2, x_vec(:).', 0, z_vec(:).', 0);
    w_out = w_out(:);
end

function gain = compute_gain_local(reference, candidate)
    denom = dot(candidate, candidate);
    if denom <= eps
        gain = 0;
    else
        gain = dot(reference, candidate) / denom;
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
    phi = zeros(nt, nz, nx, ny);
    w = zeros(nt, nz, nx, ny);
    u = zeros(nt, nz, nx, ny);
    v = zeros(nt, nz, nx, ny);
    uz = zeros(nt, nz, nx, ny);
    vz = zeros(nt, nz, nx, ny);
    wz = zeros(nt, nz, nx, ny);
    t = (0:nt - 1) * dt * tstride;

    it = 0;
    for it_read = 1:nt - 1
        tmp_eta = fread(fid, nx * ny, 'double');
        if numel(tmp_eta) < nx * ny, it = it_read - 1; break; end
        eta(it_read, :) = tmp_eta; fread(fid, 2, int_nbit);

        tmp_skip = fread(fid, nx * ny, 'double'); if numel(tmp_skip) < nx * ny, it = it_read - 1; break; end %#ok<NASGU>
        fread(fid, 2, int_nbit);
        tmp_skip = fread(fid, nx * ny, 'double'); if numel(tmp_skip) < nx * ny, it = it_read - 1; break; end %#ok<NASGU>
        fread(fid, 2, int_nbit);

        tmp_phi = fread(fid, nx * ny * nz, 'double'); if numel(tmp_phi) < nx * ny * nz, it = it_read - 1; break; end
        phi(it_read, :) = tmp_phi; fread(fid, 2, int_nbit);

        tmp_u = fread(fid, nx * ny * nz, 'double'); if numel(tmp_u) < nx * ny * nz, it = it_read - 1; break; end
        u(it_read, :) = tmp_u; fread(fid, 2, int_nbit);

        tmp_v = fread(fid, nx * ny * nz, 'double'); if numel(tmp_v) < nx * ny * nz, it = it_read - 1; break; end
        v(it_read, :) = tmp_v; fread(fid, 2, int_nbit);

        tmp_w = fread(fid, nx * ny * nz, 'double'); if numel(tmp_w) < nx * ny * nz, it = it_read - 1; break; end
        w(it_read, :) = tmp_w; fread(fid, 2, int_nbit);

        tmp_wz = fread(fid, nx * ny * nz, 'double'); if numel(tmp_wz) < nx * ny * nz, it = it_read - 1; break; end
        wz(it_read, :) = tmp_wz; fread(fid, 2, int_nbit);

        tmp_uz = fread(fid, nx * ny * nz, 'double'); if numel(tmp_uz) < nx * ny * nz, it = it_read - 1; break; end
        uz(it_read, :) = tmp_uz; fread(fid, 2, int_nbit);

        tmp_vz = fread(fid, nx * ny * nz, 'double'); if numel(tmp_vz) < nx * ny * nz, it = it_read - 1; break; end
        vz(it_read, :) = tmp_vz; fread(fid, 2, int_nbit);

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

    etat_m = 0; etatt_m = 0; phit_m = 0; p_m = 0; ut_m = 0;
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

function idx = resolve_time_index_local(CFG, it_vec, nt)
    if ~isempty(CFG.time_index)
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
    if ~isequal(size(phase_data.eta), size(ref.eta)) || ~isequal(size(phase_data.w), size(ref.w))
        error('Phase files are not compatible.');
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
