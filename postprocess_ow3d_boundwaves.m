% postprocess_ow3d_boundwaves.m
% Reconstruct first-, second-, and third-order OW3D wave components from
% four phase-shifted simulations using the same Hilbert/four-phase workflow
% as the existing extract_eta33_from_OW3D scripts.

clc;
clear;
close all;

CFG = struct();

% -------------------- User configuration --------------------
CFG.data_root = fullfile(pwd, 'directional initial condition', 'test_generator');
CFG.folder_pattern = 'T_init-10_Tp_kd1.0_spread_15_heading_0_Akp_012_alpha_8.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.time_step = 300; % EP index, e.g. 300 -> EP_00300.bin
CFG.lambda = 225;
CFG.section_mode = 'centerline'; % 'centerline', 'index', or 'y_value'
CFG.section_index = [];
CFG.section_y_value = 0.0;
CFG.save_mat = false;
CFG.output_dir = fullfile(pwd, 'processed_boundwaves');

% -------------------- Load four OW3D phase snapshots --------------------
eta_phases = [];
phi_phases = [];

for idx = 1:numel(CFG.phi_shifts_deg)
    case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
    if ~isfolder(case_folder)
        error('Missing phase folder: %s', case_folder);
    end

    bin_name = sprintf('EP_%05d.bin', CFG.time_step);
    bin_path = fullfile(case_folder, bin_name);
    if ~isfile(bin_path)
        error('Missing OW3D output: %s', bin_path);
    end

    [X, Y, eta_tmp, phi_tmp] = read_ow3d_snapshot(bin_path); %#ok<NASGU>
    eta_phases(idx,:,:) = eta_tmp;
    phi_phases(idx,:,:) = phi_tmp;
    fprintf('Loaded %s\n', bin_path);
end

% -------------------- Four-phase harmonic reconstruction -----------------
four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

[eta_harmonics, phi_harmonics] = reconstruct_harmonics(eta_phases, phi_phases, four_phase_coef);

eta1 = squeeze(eta_harmonics(1,:,:));
eta2 = squeeze(eta_harmonics(2,:,:));
eta3 = squeeze(eta_harmonics(3,:,:));
phi1 = squeeze(phi_harmonics(1,:,:));
phi2 = squeeze(phi_harmonics(2,:,:));
phi3 = squeeze(phi_harmonics(3,:,:));

kp = 2 * pi / CFG.lambda;
x_vec = X(:,1);
y_vec = Y(1,:);

% Optional harmonic cleanup in k-space for cleaner plotting/comparison.
eta1 = frequency_filtering_2d_local(eta1, x_vec, y_vec, kp, 1);
eta2 = frequency_filtering_2d_local(eta2, x_vec, y_vec, kp, 2);
eta3 = frequency_filtering_2d_local(eta3, x_vec, y_vec, kp, 3);
phi1 = frequency_filtering_2d_local(phi1, x_vec, y_vec, kp, 1);
phi2 = frequency_filtering_2d_local(phi2, x_vec, y_vec, kp, 2);
phi3 = frequency_filtering_2d_local(phi3, x_vec, y_vec, kp, 3);

fprintf('eta1 range: [%.4f, %.4f] m\n', min(eta1(:)), max(eta1(:)));
fprintf('eta2 range: [%.4f, %.4f] m\n', min(eta2(:)), max(eta2(:)));
fprintf('eta3 range: [%.4f, %.4f] m\n', min(eta3(:)), max(eta3(:)));

% -------------------- Save processed snapshot ----------------------------
if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

meta = struct();
meta.data_root = CFG.data_root;
meta.folder_pattern = CFG.folder_pattern;
meta.time_step = CFG.time_step;
meta.phi_shifts_deg = CFG.phi_shifts_deg;
meta.lambda = CFG.lambda;

if CFG.save_mat
    save(fullfile(CFG.output_dir, sprintf('OW3D_boundwaves_t%05d.mat', CFG.time_step)), ...
        'X', 'Y', 'eta1', 'eta2', 'eta3', 'phi1', 'phi2', 'phi3', ...
        'eta_harmonics', 'phi_harmonics', 'meta');
end

% -------------------- Visualization -------------------------------------
x_plot = X(:,1) - 0.5 * (X(1,1) + X(end,1));
y_plot = Y(1,:) - 0.5 * (Y(1,1) + Y(1,end));
section_idx = resolve_section_index(CFG, y_plot, size(eta1, 2));
eta_nonlinear = eta1 + eta2 + eta3;
envelope_peak_by_y = peak_envelope_along_x(eta_nonlinear);
[off_section_idx, off_section_peak] = resolve_half_height_off_section(y_plot, envelope_peak_by_y, section_idx);
fields_to_plot = {eta1, eta2, eta3, eta_nonlinear};
titles_center = { ...
    sprintf('(a) First-order elevation, centerline ($y = %.2f$ m)', y_plot(section_idx)), ...
    sprintf('(b) Second-order elevation, centerline ($y = %.2f$ m)', y_plot(section_idx)), ...
    sprintf('(c) Third-order elevation, centerline ($y = %.2f$ m)', y_plot(section_idx)), ...
    sprintf('(d) Nonlinear elevation, centerline ($y = %.2f$ m)', y_plot(section_idx))};
titles_off = { ...
    sprintf('(a) First-order elevation, off-centerline ($y = %.2f$ m)', y_plot(off_section_idx)), ...
    sprintf('(b) Second-order elevation, off-centerline ($y = %.2f$ m)', y_plot(off_section_idx)), ...
    sprintf('(c) Third-order elevation, off-centerline ($y = %.2f$ m)', y_plot(off_section_idx)), ...
    sprintf('(d) Nonlinear elevation, off-centerline ($y = %.2f$ m)', y_plot(off_section_idx))};
line_colors = [0.10 0.10 0.10; 0.80 0.26 0.18; 0.12 0.39 0.71; 0.55 0.16 0.51];
y_limits = compute_shared_ylimits(fields_to_plot, [section_idx, off_section_idx]);

fig_center = create_publishable_figure([120 100 1450 880]);
tile_center = tiledlayout(fig_center, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile_center, sprintf('OW3D bound-wave decomposition at the centerline ($t_{\\mathrm{EP}} = %d$)', CFG.time_step), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:4
    ax = nexttile(tile_center);
    draw_wave_panel(ax, x_plot, fields_to_plot{i}(:, section_idx), line_colors(i,:), titles_center{i}, y_limits(i,:));
end

annotation(fig_center, 'textbox', [0.13 0.01 0.8 0.04], ...
    'String', sprintf('Centerline defined by $y = %.2f$ m.', y_plot(section_idx)), ...
    'Interpreter', 'latex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'FontName', 'Times New Roman', 'FontSize', 11);

fig_off = create_publishable_figure([150 130 1450 880]);
tile_off = tiledlayout(fig_off, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
title(tile_off, sprintf('OW3D bound-wave decomposition at the half-envelope off-centerline ($t_{\\mathrm{EP}} = %d$)', CFG.time_step), ...
    'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:4
    ax = nexttile(tile_off);
    draw_wave_panel(ax, x_plot, fields_to_plot{i}(:, off_section_idx), line_colors(i,:), titles_off{i}, y_limits(i,:));
end

annotation(fig_off, 'textbox', [0.13 0.01 0.82 0.04], ...
    'String', sprintf('Off-centerline selected from the nonlinear envelope criterion: $y = %.2f$ m, peak envelope = %.3f m.', ...
    y_plot(off_section_idx), off_section_peak), ...
    'Interpreter', 'latex', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
    'FontName', 'Times New Roman', 'FontSize', 11);

exportgraphics(fig_center, fullfile(CFG.output_dir, sprintf('OW3D_boundwaves_centerline_t%05d.png', CFG.time_step)), 'Resolution', 300);
exportgraphics(fig_off, fullfile(CFG.output_dir, sprintf('OW3D_boundwaves_offcenterline_t%05d.png', CFG.time_step)), 'Resolution', 300);

disp('OW3D bound-wave postprocessing complete.');

% -------------------- Local helper functions -----------------------------
function [X, Y, eta_field, phi_field] = read_ow3d_snapshot(bin_path)
    [X, Y, eta_field, phi_field] = ReadBinFile_local(bin_path);
end

function [X, Y, E, P] = ReadBinFile_local(filename)
    byteorder = 'ieee-le'; % IEEE Little-Endian format
    fid = fopen(filename, 'r', byteorder);
    if fid < 0
        error('Unable to open OW3D bin file: %s', filename);
    end
    cleanup = onCleanup(@() fclose(fid));

    fread(fid, 1, 'int32');
    Nx = fread(fid, 1, 'int32');
    Ny = fread(fid, 1, 'int32');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    X = fread(fid, [Nx Ny], 'float64');
    Y = fread(fid, [Nx Ny], 'float64');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    E = fread(fid, [Nx Ny], 'float64');
    P = fread(fid, [Nx Ny], 'float64');
end

function [eta_harmonics, phi_harmonics] = reconstruct_harmonics(eta_phases, phi_phases, coef)
    eta_hilbert = hilbert2d_all(eta_phases);
    phi_hilbert = hilbert2d_all(phi_phases);

    all_eta = cat(1, real(eta_phases), -imag(eta_hilbert));
    all_phi = cat(1, real(phi_phases), -imag(phi_hilbert));

    eta_harmonics = zeros(4, size(eta_phases,2), size(eta_phases,3));
    phi_harmonics = zeros(4, size(phi_phases,2), size(phi_phases,3));

    for n = 1:4
        weights = reshape(coef(n,:), [8, 1, 1]);
        eta_harmonics(n,:,:) = sum(all_eta .* weights, 1);
        phi_harmonics(n,:,:) = sum(all_phi .* weights, 1);
    end
end

function X_hilbert = hilbert2d_all(X)
    [P, M, N] = size(X);
    X_hilbert = zeros(size(X));

    for p = 1:P
        FX = fft2(squeeze(X(p,:,:)));
        H = zeros(M, N);
        H(1:floor(M/2), :) = 1;
        H(ceil(M/2):end, :) = -1;
        X_hilbert(p,:,:) = ifft2(FX .* H);
    end
end

function field_out = frequency_filtering_2d_local(field_in, x_vec, y_vec, kp, n)
    x_vec = x_vec(:);
    y_vec = y_vec(:);
    [Nx, Ny] = size(field_in);

    if length(x_vec) ~= Nx || length(y_vec) ~= Ny
        error('Dimension mismatch between field and coordinate vectors.');
    end

    dx = x_vec(2) - x_vec(1);
    dy = y_vec(2) - y_vec(1);
    dkx = 2 * pi / (Nx * dx);
    dky = 2 * pi / (Ny * dy);

    kx = [0:ceil(Nx/2)-1, -floor(Nx/2):-1]' * dkx;
    ky = [0:ceil(Ny/2)-1, -floor(Ny/2):-1] * dky;
    [KX, KY] = ndgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);

    sigma = 0.5 * kp;
    k_target = n * kp;
    mask = exp(-((K - k_target).^2) / (2 * sigma^2));

    field_fft = fft2(field_in);
    field_out = ifft2(field_fft .* mask);
    if isreal(field_in)
        field_out = real(field_out);
    end
end

function section_idx = resolve_section_index(CFG, y_vec, ny)
    switch lower(CFG.section_mode)
        case 'centerline'
            [~, section_idx] = min(abs(y_vec));
        case 'index'
            section_idx = min(max(1, CFG.section_index), ny);
        case 'y_value'
            [~, section_idx] = min(abs(y_vec - CFG.section_y_value));
        otherwise
            error('Unsupported section_mode: %s', CFG.section_mode);
    end
end

function peak_env = peak_envelope_along_x(field_in)
    [nx, ny] = size(field_in);
    peak_env = zeros(1, ny);

    for j = 1:ny
        signal = field_in(:, j);
        analytic_signal = hilbert(signal);
        peak_env(j) = max(abs(analytic_signal(1:nx)));
    end
end

function [section_idx, section_peak] = resolve_half_height_off_section(y_vec, peak_env, center_idx)
    target_peak = 0.5 * max(peak_env);
    candidate_mask = y_vec > y_vec(center_idx);

    if ~any(candidate_mask)
        candidate_mask = true(size(y_vec));
        candidate_mask(center_idx) = false;
    end

    candidate_idx = find(candidate_mask);
    [~, best_local] = min(abs(peak_env(candidate_idx) - target_peak));
    section_idx = candidate_idx(best_local);
    section_peak = peak_env(section_idx);
end

function fig = create_publishable_figure(fig_position)
    fig = figure('Color', 'w', 'Position', fig_position, 'Renderer', 'painters');
end

function draw_wave_panel(ax, x_plot, y_plot, line_color, panel_title, y_limits)
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
    ylabel(ax, '$\eta$ (m)', 'Interpreter', 'latex', 'FontSize', 13);
    title(ax, panel_title, 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'normal');
end

function y_limits = compute_shared_ylimits(fields_to_plot, section_indices)
    n_fields = numel(fields_to_plot);
    y_limits = zeros(n_fields, 2);

    for i = 1:n_fields
        values = fields_to_plot{i}(:, section_indices);
        y_abs_max = max(abs(values(:)));
        if y_abs_max == 0
            y_abs_max = 1;
        end
        padding = 0.08 * y_abs_max;
        y_limits(i,:) = [-y_abs_max - padding, y_abs_max + padding];
    end
end
