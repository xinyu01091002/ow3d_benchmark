% postprocess_ow3d_boundkinematics.m
% Reconstruct first-, second-, and third-order OW3D bound-kinematic
% components from four phase-shifted simulations using the same
% Hilbert/four-phase workflow as the surface postprocessor.

clc;
clear;
close all;

CFG = struct();

% -------------------- User configuration --------------------
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.kinematics_file_id = 1; % Kinematics01.bin
CFG.time_index = []; % [] -> use default near-final frame. Positive -> index from start. Negative -> index from end.
CFG.default_time_index_from_end = 160; % Used only when time_index = [].
CFG.lambda = 225;
CFG.variables_to_process = {'u', 'w', 'phi', 'p'};
CFG.apply_x_filter = true;
CFG.sigma_mode = 'surface'; % 'surface', 'index', or 'value'
CFG.sigma_index = [];
CFG.sigma_value = 0.0;
CFG.save_mat = false;
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');

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

    phase_data = read_ow3d_kinematics_snapshot(kin_path);
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

fprintf('Using kinematics time index %d of %d (t = %.6f s)\n', ...
    selected_time_index, numel(t_vec), t_selected);

eta_phases = zeros(numel(CFG.phi_shifts_deg), size(ref.eta, 2));
vars_phases = initialize_variable_phase_storage(CFG.variables_to_process, numel(CFG.phi_shifts_deg), ref);

for idx = 1:numel(CFG.phi_shifts_deg)
    phase_data = data_by_phase{idx};
    eta_phases(idx, :) = squeeze(phase_data.eta(selected_time_index, :, 1));

    for v_idx = 1:numel(CFG.variables_to_process)
        var_name = CFG.variables_to_process{v_idx};
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

for v_idx = 1:numel(CFG.variables_to_process)
    var_name = CFG.variables_to_process{v_idx};
    var_harmonics.(var_name) = reconstruct_harmonics_xz(vars_phases.(var_name), four_phase_coef);
end

% -------------------- Optional x-direction harmonic cleanup --------------
kp = 2 * pi / CFG.lambda;
if CFG.apply_x_filter
    eta_harmonics = filter_harmonics_x_only(eta_harmonics, x_vec, kp);

    for v_idx = 1:numel(CFG.variables_to_process)
        var_name = CFG.variables_to_process{v_idx};
        var_harmonics.(var_name) = filter_harmonics_x_only(var_harmonics.(var_name), x_vec, kp);
    end
end

% -------------------- Save processed snapshot ----------------------------
if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

meta = struct();
meta.data_root = CFG.data_root;
meta.folder_pattern = CFG.folder_pattern;
meta.phi_shifts_deg = CFG.phi_shifts_deg;
meta.kinematics_file_id = CFG.kinematics_file_id;
meta.time_index = selected_time_index;
meta.time_value = t_selected;
meta.lambda = CFG.lambda;
meta.sigma = sigma_vec;
meta.x = x_vec;
meta.y = y_vec;

if CFG.save_mat
    save(fullfile(CFG.output_dir, sprintf('OW3D_boundkinematics_tidx_%04d.mat', selected_time_index)), ...
        'eta_harmonics', 'var_harmonics', 'meta');
end

% -------------------- Visualization -------------------------------------
x_plot = x_vec - 0.5 * (x_vec(1) + x_vec(end));
sigma_idx = resolve_sigma_index(CFG, sigma_vec);
sigma_value = sigma_vec(sigma_idx);
line_colors = [0.10 0.10 0.10; 0.80 0.26 0.18; 0.12 0.39 0.71; 0.55 0.16 0.51];

for v_idx = 1:numel(CFG.variables_to_process)
    var_name = CFG.variables_to_process{v_idx};
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
        draw_profile_panel(ax, x_plot, fields_to_plot{panel_idx}, line_colors(panel_idx, :), titles{panel_idx}, y_limits(panel_idx, :), variable_axis_label(var_name));
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

disp('OW3D bound-kinematics postprocessing complete.');

% -------------------- Local helper functions -----------------------------
function data = read_ow3d_kinematics_snapshot(kin_path)
    [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, sigma, t] = ...
        read_kinematics_file_local(kin_path); %#ok<ASGLU>

    data = struct();
    data.it = it;
    data.eta = eta;
    data.etat_m = etat_m;
    data.etatt_m = etatt_m;
    data.phi = phi;
    data.phit_m = phit_m;
    data.p = p_m;
    data.ut = ut_m;
    data.u = u;
    data.v = v;
    data.w = w;
    data.uz = uz;
    data.vz = vz;
    data.wz = wz;
    data.x = x;
    data.y = y;
    data.sigma = sigma;
    data.t = t;
end

function [it, eta, etat_m, etatt_m, phi, phit_m, p_m, ut_m, u, v, w, uz, vz, wz, x, y, sigma, t] = read_kinematics_file_local(file_path)
    nbits = 32;
    compute_pressure = true;

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

    if compute_pressure
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

                    phit(:, j, ip) = dt_matrix * phi_col - w_col .* sigma(j) .* etat(:, ip);
                    p(:, j, ip) = -(phit(:, j, ip) + 0.5 * (u_col.^2 + w_col.^2));
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

function draw_profile_panel(ax, x_plot, y_plot, line_color, panel_title, y_limits, y_label)
    plot(ax, x_plot, y_plot, 'Color', line_color, 'LineWidth', 1.8);
    hold(ax, 'on');
    yline(ax, 0, '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
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
    xlim(ax, [x_plot(1), x_plot(end)]);
    ylim(ax, y_limits);
    xlabel(ax, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
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
        case 'p'
            label = 'dynamic pressure surrogate';
        case 'ut'
            label = 'horizontal acceleration';
        case 'phit'
            label = 'potential time derivative';
        otherwise
            label = var_name;
    end
end

function label = variable_axis_label(var_name)
    switch lower(var_name)
        case {'u', 'v', 'w', 'ut'}
            label = sprintf('$%s$', var_name);
        case 'phi'
            label = '$\phi$';
        case 'p'
            label = '$p$';
        case 'phit'
            label = '$\phi_t$';
        otherwise
            label = ['$', var_name, '$'];
    end
end
