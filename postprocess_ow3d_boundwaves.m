% postprocess_ow3d_boundwaves.m
% Reconstruct first-, second-, and third-order OW3D wave components from
% four phase-shifted simulations using the same Hilbert/four-phase workflow
% as the existing extract_eta33_from_OW3D scripts.

clc;
clear;
close all;

CFG = struct();

% -------------------- User configuration --------------------
CFG.data_root = fullfile(pwd, 'directional initial condition', 'error_wave_separation');
CFG.folder_pattern = 'T_init-40_Tend5_Tp_kd1.0_spread_15_heading_0_Akp_012_alpha_8.0_phi_%d';
CFG.phi_shifts_deg = 0:90:270;
CFG.lambda = 225;
CFG.section_mode = 'centerline'; % 'centerline', 'index', or 'y_value'
CFG.section_index = [];
CFG.section_y_value = 0.0;
CFG.save_mat = false;
CFG.output_dir = fullfile(pwd, 'processed_boundwaves');
CFG.export_initial_snapshot = true;
CFG.export_peak_envelope_snapshot = true;
CFG.initial_time_step = []; % Empty -> use the first available EP file in the phase-0 folder.
CFG.case_tag_override = ''; % Optional short tag for filenames. Empty -> auto derived from folder_pattern.
CFG.ignore_time_steps = [99999]; % Ignore known placeholder/sentinel EP files.
CFG.ylim_mode = 'per_panel'; % 'per_panel', 'per_order', or 'global'
CFG.ylim_padding_fraction = 0.10;
CFG.ylim_min_abs = 1e-3; % Minimum half-range so very small signals remain visible.

if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

phase0_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(1)));
if ~isfolder(phase0_folder)
    error('Missing phase-0 folder: %s', phase0_folder);
end

case_tag = resolve_case_tag(CFG);
target_steps = determine_target_time_steps(CFG, phase0_folder);

for target_idx = 1:numel(target_steps)
    target = target_steps(target_idx);
    fprintf('Processing %s snapshot at EP_%05d.bin\n', target.label, target.time_step);

    [X, Y, eta1, eta2, eta3, phi1, phi2, phi3, eta_harmonics, phi_harmonics] = ...
        process_single_time_step(CFG, target.time_step);

    fprintf('  eta1 range: [%.4f, %.4f] m\n', min(eta1(:)), max(eta1(:)));
    fprintf('  eta2 range: [%.4f, %.4f] m\n', min(eta2(:)), max(eta2(:)));
    fprintf('  eta3 range: [%.4f, %.4f] m\n', min(eta3(:)), max(eta3(:)));

    meta = struct();
    meta.data_root = CFG.data_root;
    meta.folder_pattern = CFG.folder_pattern;
    meta.time_step = target.time_step;
    meta.time_label = target.label;
    meta.phi_shifts_deg = CFG.phi_shifts_deg;
    meta.lambda = CFG.lambda;
    meta.case_tag = case_tag;

    if CFG.save_mat
        save(fullfile(CFG.output_dir, sprintf('OW3D_boundwaves_%s_%s_t%05d.mat', ...
            case_tag, target.label, target.time_step)), ...
            'X', 'Y', 'eta1', 'eta2', 'eta3', 'phi1', 'phi2', 'phi3', ...
            'eta_harmonics', 'phi_harmonics', 'meta');
    end

    export_boundwave_figures(CFG, X, Y, eta1, eta2, eta3, target, case_tag);
end

disp('OW3D bound-wave postprocessing complete.');

% -------------------- Local helper functions -----------------------------
function [X, Y, eta1, eta2, eta3, phi1, phi2, phi3, eta_harmonics, phi_harmonics] = process_single_time_step(CFG, time_step)
    eta_phases = [];
    phi_phases = [];

    for idx = 1:numel(CFG.phi_shifts_deg)
        case_folder = fullfile(CFG.data_root, sprintf(CFG.folder_pattern, CFG.phi_shifts_deg(idx)));
        if ~isfolder(case_folder)
            error('Missing phase folder: %s', case_folder);
        end

        bin_name = sprintf('EP_%05d.bin', time_step);
        bin_path = fullfile(case_folder, bin_name);
        if ~isfile(bin_path)
            error('Missing OW3D output: %s', bin_path);
        end

        [X, Y, eta_tmp, phi_tmp] = read_ow3d_snapshot(bin_path); %#ok<AGROW>
        eta_phases(idx,:,:) = eta_tmp; %#ok<AGROW>
        phi_phases(idx,:,:) = phi_tmp; %#ok<AGROW>
    end

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

    eta1 = frequency_filtering_2d_local(eta1, x_vec, y_vec, kp, 1);
    eta2 = frequency_filtering_2d_local(eta2, x_vec, y_vec, kp, 2);
    eta3 = frequency_filtering_2d_local(eta3, x_vec, y_vec, kp, 3);
    phi1 = frequency_filtering_2d_local(phi1, x_vec, y_vec, kp, 1);
    phi2 = frequency_filtering_2d_local(phi2, x_vec, y_vec, kp, 2);
    phi3 = frequency_filtering_2d_local(phi3, x_vec, y_vec, kp, 3);
end

function export_boundwave_figures(CFG, X, Y, eta1, eta2, eta3, target, case_tag)
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
    center_profiles = {
        fields_to_plot{1}(:, section_idx), ...
        fields_to_plot{2}(:, section_idx), ...
        fields_to_plot{3}(:, section_idx), ...
        fields_to_plot{4}(:, section_idx)};
    off_profiles = {
        fields_to_plot{1}(:, off_section_idx), ...
        fields_to_plot{2}(:, off_section_idx), ...
        fields_to_plot{3}(:, off_section_idx), ...
        fields_to_plot{4}(:, off_section_idx)};
    y_limits_center = compute_dynamic_ylimits(CFG, center_profiles, off_profiles);
    y_limits_off = compute_dynamic_ylimits(CFG, off_profiles, center_profiles);

    fig_center = create_publishable_figure([120 80 1500 980]);
    tile_center = tiledlayout(fig_center, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    configure_tiled_layout(tile_center);
    title(tile_center, sprintf('OW3D bound-wave decomposition at the centerline (%s, $t_{\\mathrm{EP}} = %d$)', ...
        format_time_label(target.label), target.time_step), ...
        'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

    for i = 1:4
        ax = nexttile(tile_center);
        draw_wave_panel(ax, x_plot, center_profiles{i}, line_colors(i,:), titles_center{i}, y_limits_center(i,:));
    end

    add_footer_note(fig_center, sprintf('Centerline defined by $y = %.2f$ m.', y_plot(section_idx)));

    fig_off = create_publishable_figure([150 110 1500 980]);
    tile_off = tiledlayout(fig_off, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    configure_tiled_layout(tile_off);
    title(tile_off, sprintf('OW3D bound-wave decomposition at the half-envelope off-centerline (%s, $t_{\\mathrm{EP}} = %d$)', ...
        format_time_label(target.label), target.time_step), ...
        'Interpreter', 'latex', 'FontSize', 15, 'FontWeight', 'bold');

    for i = 1:4
        ax = nexttile(tile_off);
        draw_wave_panel(ax, x_plot, off_profiles{i}, line_colors(i,:), titles_off{i}, y_limits_off(i,:));
    end

    add_footer_note(fig_off, sprintf(['Off-centerline selected from the nonlinear envelope criterion: ' ...
        '$y = %.2f$ m, peak envelope = %.3f m.'], y_plot(off_section_idx), off_section_peak));

    exportgraphics(fig_center, fullfile(CFG.output_dir, build_export_name(case_tag, target.label, target.time_step, 'centerline')), 'Resolution', 300);
    exportgraphics(fig_off, fullfile(CFG.output_dir, build_export_name(case_tag, target.label, target.time_step, 'offcenterline')), 'Resolution', 300);
end

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

function configure_tiled_layout(tile)
    tile.Units = 'normalized';
    tile.Position = [0.08 0.12 0.88 0.80];
end

function add_footer_note(fig, note_text)
    annotation(fig, 'textbox', [0.10 0.025 0.84 0.05], ...
        'String', note_text, ...
        'Interpreter', 'latex', ...
        'EdgeColor', 'none', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 11, ...
        'FitBoxToText', 'off');
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
    title(ax, panel_title, 'Interpreter', 'latex', 'FontSize', 12, 'FontWeight', 'normal');
end

function y_limits = compute_dynamic_ylimits(CFG, primary_profiles, companion_profiles)
    n_fields = numel(primary_profiles);
    y_limits = zeros(n_fields, 2);

    for i = 1:n_fields
        primary_abs_max = max(abs(primary_profiles{i}(:)));

        switch lower(CFG.ylim_mode)
            case 'per_panel'
                ref_abs_max = primary_abs_max;
            case 'per_order'
                companion_abs_max = max(abs(companion_profiles{i}(:)));
                ref_abs_max = max(primary_abs_max, companion_abs_max);
            case 'global'
                ref_abs_max = primary_abs_max;
                for j = 1:n_fields
                    ref_abs_max = max(ref_abs_max, max(abs(primary_profiles{j}(:))));
                    ref_abs_max = max(ref_abs_max, max(abs(companion_profiles{j}(:))));
                end
            otherwise
                error('Unsupported ylim_mode: %s', CFG.ylim_mode);
        end

        ref_abs_max = max(ref_abs_max, CFG.ylim_min_abs);
        padding = CFG.ylim_padding_fraction * ref_abs_max;
        y_limits(i,:) = [-ref_abs_max - padding, ref_abs_max + padding];
    end
end

function target_steps = determine_target_time_steps(CFG, phase0_folder)
    available_steps = list_available_time_steps(CFG, phase0_folder);
    if isempty(available_steps)
        error('No EP_*.bin files found in %s', phase0_folder);
    end

    target_steps = struct('label', {}, 'time_step', {});

    if CFG.export_initial_snapshot
        initial_step = resolve_initial_time_step(CFG.initial_time_step, available_steps);
        target_steps(end + 1).label = 'initial'; %#ok<AGROW>
        target_steps(end).time_step = initial_step;
    end

    if CFG.export_peak_envelope_snapshot
        peak_step = find_peak_envelope_time_step(phase0_folder, available_steps);
        if isempty(target_steps) || peak_step ~= target_steps(end).time_step
            target_steps(end + 1).label = 'peakEnvelope'; %#ok<AGROW>
            target_steps(end).time_step = peak_step;
        else
            target_steps(end).label = 'initial_peakEnvelope';
        end
    end

    if isempty(target_steps)
        error('No export target selected. Enable initial and/or peak-envelope snapshot export.');
    end
end

function available_steps = list_available_time_steps(CFG, case_folder)
    files = dir(fullfile(case_folder, 'EP_*.bin'));
    raw_steps = zeros(1, numel(files));

    for i = 1:numel(files)
        token = regexp(files(i).name, '^EP_(\d+)\.bin$', 'tokens', 'once');
        if ~isempty(token)
            raw_steps(i) = str2double(token{1});
        end
    end

    available_steps = sort(unique(raw_steps(raw_steps >= 0)));
    available_steps = filter_suspicious_time_steps(available_steps, CFG.ignore_time_steps);
end

function initial_step = resolve_initial_time_step(cfg_initial_time_step, available_steps)
    if isempty(cfg_initial_time_step)
        initial_step = available_steps(1);
    else
        [~, idx] = min(abs(available_steps - cfg_initial_time_step));
        initial_step = available_steps(idx);
    end
end

function peak_step = find_peak_envelope_time_step(case_folder, available_steps)
    peak_metric = -inf;
    peak_step = available_steps(1);

    for i = 1:numel(available_steps)
        bin_path = fullfile(case_folder, sprintf('EP_%05d.bin', available_steps(i)));
        [~, ~, eta_field, ~] = read_ow3d_snapshot(bin_path);
        envelope_peak_by_y = peak_envelope_along_x(eta_field);
        current_metric = max(envelope_peak_by_y);

        if current_metric > peak_metric
            peak_metric = current_metric;
            peak_step = available_steps(i);
        end
    end
end

function case_tag = resolve_case_tag(CFG)
    if ~isempty(CFG.case_tag_override)
        case_tag = sanitize_filename_token(CFG.case_tag_override);
        return;
    end

    tag_source = strrep(CFG.folder_pattern, '_phi_%d', '');
    tag_source = strrep(tag_source, '%d', '');
    case_tag = sanitize_filename_token(tag_source);
end

function file_name = build_export_name(case_tag, time_label, time_step, section_tag)
    file_name = sprintf('OW3D_boundwaves_%s_%s_t%05d_%s.png', ...
        sanitize_filename_token(case_tag), ...
        sanitize_filename_token(time_label), ...
        time_step, ...
        sanitize_filename_token(section_tag));
end

function label_out = format_time_label(label_in)
    switch label_in
        case 'initial'
            label_out = 'initial snapshot';
        case 'peakEnvelope'
            label_out = 'peak-envelope snapshot';
        case 'initial_peakEnvelope'
            label_out = 'initial / peak-envelope snapshot';
        otherwise
            label_out = strrep(label_in, '_', '\_');
    end
end

function token = sanitize_filename_token(token_in)
    token = regexprep(token_in, '[^\w\-]+', '_');
    token = regexprep(token, '_+', '_');
    token = regexprep(token, '^_|_$', '');
    if isempty(token)
        token = 'case';
    end
end

function filtered_steps = filter_suspicious_time_steps(available_steps, ignore_time_steps)
    filtered_steps = available_steps;

    if isempty(filtered_steps)
        return;
    end

    filtered_steps = filtered_steps(~ismember(filtered_steps, ignore_time_steps));

    if numel(filtered_steps) <= 2
        return;
    end

    step_diff = diff(filtered_steps);
    positive_diff = step_diff(step_diff > 0);
    if isempty(positive_diff)
        return;
    end

    typical_diff = median(positive_diff);
    if typical_diff <= 0
        return;
    end

    last_gap = filtered_steps(end) - filtered_steps(end - 1);
    if last_gap > 10 * typical_diff
        filtered_steps(end) = [];
    end
end
