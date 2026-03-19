% test_kinematics_phase_sanity.m
% Quick sanity checks for OW3D kinematics files before applying any
% four-phase separation. This script compares the raw phase-shifted
% kinematics directly at a fixed sigma level, while also plotting the
% surface elevation as a reference.

clc;
clear;
close all;

CFG = struct();

% -------------------- User configuration --------------------
CFG.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check');
CFG.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
CFG.phases_deg = [0, 90, 180,270]; % Can be three or four available phase cases.
CFG.kinematics_file_id = 1; % Kinematics01.bin
CFG.time_index = 12; % Avoid the first stored frame if it is a zeroed startup frame.
CFG.sigma_mode = 'surface'; % 'surface', 'index', or 'value'
CFG.sigma_index = [];
CFG.sigma_value = 0.0;
CFG.variables_to_plot = {'eta', 'phi', 'u', 'w'};
CFG.output_dir = fullfile(pwd, 'processed_boundkinematics');
CFG.save_png = true;

if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

phase_data = cell(1, numel(CFG.phases_deg));
phase_labels = strings(1, numel(CFG.phases_deg));

for idx = 1:numel(CFG.phases_deg)
    phase_deg = CFG.phases_deg(idx);
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, phase_deg));
    if ~isfolder(case_folder)
        error('Missing phase folder: %s', case_folder);
    end

    kin_path = resolve_kinematics_path(case_folder, CFG.kinematics_file_id);
    phase_data{idx} = read_kinematics_file_local(kin_path);
    phase_labels(idx) = sprintf('\\phi = %d^\\circ', phase_deg);
    fprintf('Loaded %s\n', kin_path);
end

assert_phase_family_consistency(phase_data);

ref = phase_data{1};
time_index = min(max(1, CFG.time_index), ref.it);
t_value = ref.t(time_index);
x_vec = ref.x(:, 1);
sigma_vec = ref.sigma(:);
sigma_idx = resolve_sigma_index(CFG, sigma_vec);
sigma_value = sigma_vec(sigma_idx);

fprintf('Using time index %d (t = %.6f s)\n', time_index, t_value);
fprintf('Using sigma index %d (sigma = %.6f)\n', sigma_idx, sigma_value);

for v_idx = 1:numel(CFG.variables_to_plot)
    var_name = CFG.variables_to_plot{v_idx};
    fig = figure('Color', 'w', 'Position', [100 100 1450 880], 'Renderer', 'painters');
    tile = tiledlayout(fig, 1, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('Raw OW3D %s across phase cases at t = %.4f s, \\sigma = %.4f', ...
        variable_display_name(var_name), t_value, sigma_value), ...
        'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    ax = nexttile(tile);
    hold(ax, 'on');
    line_colors = lines(numel(CFG.phases_deg));

    for p_idx = 1:numel(CFG.phases_deg)
        profile = extract_x_profile(phase_data{p_idx}, var_name, time_index, sigma_idx);
        plot(ax, x_vec, profile, 'LineWidth', 1.6, ...
            'Color', line_colors(p_idx, :), ...
            'DisplayName', phase_labels(p_idx));
    end

    hold(ax, 'off');
    grid(ax, 'on');
    box(ax, 'on');
    ax.LineWidth = 1.0;
    ax.FontName = 'Times New Roman';
    ax.FontSize = 12;
    xlabel(ax, '$x$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
    ylabel(ax, variable_axis_label(var_name), 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, sprintf('%s profile', variable_display_name(var_name)), 'Interpreter', 'tex', ...
        'FontSize', 13, 'FontWeight', 'normal');
    legend(ax, 'Location', 'best', 'Box', 'off');

    if CFG.save_png
        exportgraphics(fig, fullfile(CFG.output_dir, sprintf('kinematics_phase_profile_%s_sigma_%03d_sanity.png', var_name, sigma_idx)), 'Resolution', 300);
    end
end

if numel(CFG.phases_deg) >= 3 && any(CFG.phases_deg == 0) && any(CFG.phases_deg == 180)
    idx0 = find(CFG.phases_deg == 0, 1, 'first');
    idx180 = find(CFG.phases_deg == 180, 1, 'first');
    x_idx = round(0.5 * (numel(x_vec) + 1));
    fprintf('\nSimple 0/180 checks at time index %d, sigma index %d, x index %d:\n', time_index, sigma_idx, x_idx);

    for v_idx = 1:numel(CFG.variables_to_plot)
        var_name = CFG.variables_to_plot{v_idx};
        v0 = extract_point_value(phase_data{idx0}, var_name, time_index, sigma_idx, x_idx);
        v180 = extract_point_value(phase_data{idx180}, var_name, time_index, sigma_idx, x_idx);
        fprintf('  %-4s : value(phi=0) = %+ .6e, value(phi=180) = %+ .6e, sum = %+ .6e\n', ...
            var_name, v0, v180, v0 + v180);
    end
end

disp('Kinematics fixed-sigma phase sanity check complete.');

% -------------------- Local helper functions -----------------------------
function data = read_kinematics_file_local(file_path)
    nbits = 32;

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
    x(:) = tmp(1:5:5 * nx * ny);
    y(:) = tmp(2:5:5 * nx * ny);

    for i = 1:nz
        sigma(i) = fread(fid, 1, 'double');
    end
    fread(fid, 2, int_nbit);

    eta = zeros(nt, nx, ny);
    phi = zeros(nt, nz, nx, ny);
    u = zeros(nt, nz, nx, ny);
    v = zeros(nt, nz, nx, ny);
    w = zeros(nt, nz, nx, ny);
    uz = zeros(nt, nz, nx, ny);
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

        fread(fid, nx * ny, 'double'); % etax
        fread(fid, 2, int_nbit);
        fread(fid, nx * ny, 'double'); % etay
        fread(fid, 2, int_nbit);

        tmp_phi = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_phi) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        phi(it_read, :) = tmp_phi;
        fread(fid, 2, int_nbit);

        tmp_u = fread(fid, nx * ny * nz, 'double');
        tmp_v = fread(fid, nx * ny * nz, 'double');
        tmp_w = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_u) < nx * ny * nz || numel(tmp_v) < nx * ny * nz || numel(tmp_w) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        u(it_read, :) = tmp_u;
        fread(fid, 2, int_nbit);
        v(it_read, :) = tmp_v;
        fread(fid, 2, int_nbit);
        w(it_read, :) = tmp_w;
        fread(fid, 2, int_nbit);

        tmp_wz = fread(fid, nx * ny * nz, 'double');
        tmp_uz = fread(fid, nx * ny * nz, 'double');
        tmp_vz = fread(fid, nx * ny * nz, 'double');
        if numel(tmp_wz) < nx * ny * nz || numel(tmp_uz) < nx * ny * nz || numel(tmp_vz) < nx * ny * nz
            it = it_read - 1;
            break;
        end
        wz(it_read, :) = tmp_wz;
        fread(fid, 2, int_nbit);
        uz(it_read, :) = tmp_uz;
        fread(fid, 2, int_nbit);
        vz(it_read, :) = tmp_vz;
        fread(fid, 2, int_nbit);

        it = it_read;
    end

    if it <= 0
        error('No complete stored kinematics time step could be read from %s', file_path);
    end

    data = struct();
    data.it = it;
    data.eta = eta(1:it, :, :);
    data.phi = phi(1:it, :, :, :);
    data.u = u(1:it, :, :, :);
    data.v = v(1:it, :, :, :);
    data.w = w(1:it, :, :, :);
    data.uz = uz(1:it, :, :, :);
    data.vz = vz(1:it, :, :, :);
    data.wz = wz(1:it, :, :, :);
    data.x = x;
    data.y = y;
    data.sigma = sigma;
    data.t = t(1:it);
end

function kin_path = resolve_kinematics_path(case_folder, file_id)
    if file_id < 10
        file_name = sprintf('Kinematics0%d.bin', file_id);
    else
        file_name = sprintf('Kinematics%d.bin', file_id);
    end

    kin_path = fullfile(case_folder, file_name);
end

function assert_phase_family_consistency(phase_data)
    ref = phase_data{1};

    for idx = 2:numel(phase_data)
        cur = phase_data{idx};
        if ~isequal(size(cur.x), size(ref.x)) || ~isequal(size(cur.y), size(ref.y))
            error('Grid mismatch between phase datasets.');
        end
        if ~isequal(size(cur.eta), size(ref.eta)) || ~isequal(size(cur.u), size(ref.u))
            error('Field size mismatch between phase datasets.');
        end
        if numel(cur.sigma) ~= numel(ref.sigma) || any(abs(cur.sigma(:) - ref.sigma(:)) > 1e-12)
            error('Sigma mismatch between phase datasets.');
        end
        if numel(cur.t) ~= numel(ref.t) || any(abs(cur.t(:) - ref.t(:)) > 1e-12)
            error('Stored time mismatch between phase datasets.');
        end
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

function profile = extract_x_profile(data, var_name, time_index, sigma_idx)
    switch lower(var_name)
        case 'eta'
            profile = squeeze(data.eta(time_index, :, 1));
        otherwise
            profile = squeeze(data.(var_name)(time_index, sigma_idx, :, 1));
    end
    profile = profile(:);
end

function value = extract_point_value(data, var_name, time_index, sigma_idx, x_idx)
    switch lower(var_name)
        case 'eta'
            value = data.eta(time_index, x_idx, 1);
        otherwise
            value = data.(var_name)(time_index, sigma_idx, x_idx, 1);
    end
end

function label = variable_display_name(var_name)
    switch lower(var_name)
        case 'eta'
            label = 'surface elevation';
        case 'u'
            label = 'horizontal velocity';
        case 'v'
            label = 'transverse velocity';
        case 'w'
            label = 'vertical velocity';
        case 'phi'
            label = 'velocity potential';
        otherwise
            label = var_name;
    end
end

function label = variable_axis_label(var_name)
    switch lower(var_name)
        case 'eta'
            label = '$\eta$';
        case {'u', 'v', 'w'}
            label = sprintf('$%s$', var_name);
        case 'phi'
            label = '$\phi$';
        otherwise
            label = ['$', var_name, '$'];
    end
end
