% compare_ow3d_mf12_w20_vertical_profile.m
% Plot vertical profiles of OW3D and MF12 w20 near the free surface.

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
CFG.top_sigma_layers = 8;
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

disp('W20 vertical profile plots complete.');

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

    coeffs2 = build_mf12_coeffs2_local(eta11_surface, x_vec, depth_value, CFG.gravity, CFG.linear_energy_keep);

    nz = numel(sigma_vec);
    sigma_start = max(1, nz - CFG.top_sigma_layers + 1);
    sigma_indices = sigma_start:nz;

    ow3d_w20_layers = zeros(numel(sigma_indices), numel(x_vec));
    z_mean_layers = zeros(numel(sigma_indices), 1);
    z_moving_layers = zeros(numel(sigma_indices), numel(x_vec));

    for i = 1:numel(sigma_indices)
        sigma_idx = sigma_indices(i);
        sigma_value = sigma_vec(sigma_idx);
        w_surface_phases = squeeze(w_phases(:, sigma_idx, :));
        w_harmonics = reconstruct_harmonics_1d_local(w_surface_phases, four_phase_coef);
        ow3d_w20_layers(i, :) = w_harmonics(4, :);
        z_mean_layers(i) = depth_value * (sigma_value - 1);
        z_moving_layers(i, :) = z_mean_layers(i) + sigma_value * eta11_surface.';
    end

    [~, x_peak_idx] = max(abs(ow3d_w20_layers(end, :)));
    x_peak = x_vec(x_peak_idx);
    eta1_peak = eta11_surface(x_peak_idx);
    z_sample = linspace(min(z_mean_layers) * 1.05, max(eta11_surface) * 1.05, 240).';
    mf12_profile = evaluate_mf12_w20_profile_local(coeffs2, x_peak, z_sample);

    ow3d_profile_mean = ow3d_w20_layers(:, x_peak_idx);
    ow3d_profile_moving_z = z_moving_layers(:, x_peak_idx);

    fig = figure('Color', 'w', 'Position', [140 100 1450 880], 'Renderer', 'painters');
    tile = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('%s: w_{20} vertical structure near the free surface (t index = %d, t = %.4f s)', ...
        strrep(case_tag, '_', '\_'), selected_time_index, t_selected), ...
        'Interpreter', 'tex', 'FontSize', 15, 'FontWeight', 'bold');

    ax1 = nexttile(tile, 1);
    hold(ax1, 'on');
    plot(ax1, mf12_profile, z_sample, 'k-', 'LineWidth', 1.9, 'DisplayName', 'MF12 continuous w_{20}(z)');
    plot(ax1, ow3d_profile_mean, z_mean_layers, 'o', 'MarkerSize', 7, 'LineWidth', 1.4, 'DisplayName', 'OW3D top-layer points at mean z');
    plot(ax1, ow3d_profile_mean, ow3d_profile_moving_z, 's', 'MarkerSize', 7, 'LineWidth', 1.4, 'DisplayName', 'OW3D top-layer points at moving z');
    hold(ax1, 'off');
    grid(ax1, 'on'); box(ax1, 'on');
    xlabel(ax1, '$w_{20}$ (m/s)', 'Interpreter', 'latex');
    ylabel(ax1, '$z$ (m)', 'Interpreter', 'latex');
    title(ax1, sprintf('Profile at x = %.3f m (surface |w_{20}| max)', x_peak), 'Interpreter', 'tex');
    legend(ax1, 'Location', 'best', 'FontSize', 9);

    ax2 = nexttile(tile, 2);
    hold(ax2, 'on');
    plot(ax2, sigma_vec(sigma_indices), ow3d_profile_mean, 'o-', 'LineWidth', 1.6, 'DisplayName', 'OW3D w_{20}');
    plot(ax2, sigma_vec(sigma_indices), interp1(z_sample, mf12_profile, z_mean_layers, 'linear', 'extrap'), 'd-', ...
        'LineWidth', 1.6, 'DisplayName', 'MF12 at z=h(\sigma-1)');
    plot(ax2, sigma_vec(sigma_indices), interp1(z_sample, mf12_profile, ow3d_profile_moving_z, 'linear', 'extrap'), 's-', ...
        'LineWidth', 1.6, 'DisplayName', 'MF12 at z=\sigma(h+\eta^{(1)})-h');
    hold(ax2, 'off');
    grid(ax2, 'on'); box(ax2, 'on');
    xlabel(ax2, '\sigma'); ylabel(ax2, '$w_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax2, 'Top-layer amplitudes at the selected x');
    legend(ax2, 'Location', 'best', 'FontSize', 9);

    ax3 = nexttile(tile, 3);
    hold(ax3, 'on');
    plot(ax3, x_vec, ow3d_w20_layers(end, :), 'k-', 'LineWidth', 1.8, 'DisplayName', sprintf('OW3D w_{20}, \\sigma=%.4f', sigma_vec(end)));
    xline(ax3, x_peak, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', 'selected x');
    hold(ax3, 'off');
    grid(ax3, 'on'); box(ax3, 'on');
    xlabel(ax3, '$x$ (m)', 'Interpreter', 'latex');
    ylabel(ax3, '$w_{20}$ (m/s)', 'Interpreter', 'latex');
    title(ax3, 'Surface-layer w_{20} used to select x');
    legend(ax3, 'Location', 'best', 'FontSize', 9);

    ax4 = nexttile(tile, 4);
    hold(ax4, 'on');
    plot(ax4, x_vec, eta11_surface, 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6, 'DisplayName', '\eta^{(1)}');
    xline(ax4, x_peak, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.2, 'DisplayName', 'selected x');
    hold(ax4, 'off');
    grid(ax4, 'on'); box(ax4, 'on');
    xlabel(ax4, '$x$ (m)', 'Interpreter', 'latex');
    ylabel(ax4, '$\eta^{(1)}$ (m)', 'Interpreter', 'latex');
    title(ax4, sprintf('Selected x and local \\eta^{(1)} = %.4f m', eta1_peak), 'Interpreter', 'tex');
    legend(ax4, 'Location', 'best', 'FontSize', 9);

    exportgraphics(fig, fullfile(output_dir, sprintf('compare_ow3d_mf12_w20_vertical_profile_tidx_%04d.png', selected_time_index)), 'Resolution', 300);

    results = struct();
    results.case_tag = case_tag;
    results.t_index = selected_time_index;
    results.t_value = t_selected;
    results.x_peak_idx = x_peak_idx;
    results.x_peak = x_peak;
    results.eta1_peak = eta1_peak;
    results.sigma_indices = sigma_indices;
    results.sigma_values = sigma_vec(sigma_indices);
    results.z_mean_layers = z_mean_layers;
    results.z_moving_layers_at_peak = ow3d_profile_moving_z;
    results.ow3d_w20_layers = ow3d_w20_layers(:, x_peak_idx);
    results.mf12_profile_z = z_sample;
    results.mf12_profile_w20 = mf12_profile;
    save(fullfile(output_dir, sprintf('compare_ow3d_mf12_w20_vertical_profile_tidx_%04d.mat', selected_time_index)), 'results', '-v7.3');
end

function coeffs2 = build_mf12_coeffs2_local(eta11, x_vec, depth, gravity, energy_keep)
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
end

function w_profile = evaluate_mf12_w20_profile_local(coeffs2, x0, z_vec)
    npts = numel(z_vec);
    w_profile = zeros(npts, 1);
    for idx = 1:npts
        [~, ~, ~, ~, w_profile(idx)] = mf12_second_subharmonic_kinematics(coeffs2, x0, 0, z_vec(idx), 0);
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
