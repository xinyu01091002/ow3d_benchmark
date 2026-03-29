% make_directional_harmonic_evolution_video.m
% Reconstruct four-phase directional OW3D surface harmonics and export
% a 2x2 video showing the evolution of:
%   (1) first-order component
%   (2) second-order superharmonic
%   (3) third-order superharmonic
%   (4) second-order subharmonic / difference-frequency component

clc;
clear;
close all;

CFG = struct();
CFG.data_root = fullfile(pwd, 'directional initial condition', 'cpp_large_spreading');
CFG.output_dir = fullfile(pwd, 'videos');
CFG.output_name = 'directional_cpp_large_spreading_harmonic_evolution.mp4';
CFG.centerline_output_name = 'directional_cpp_large_spreading_eta11_centerline_linear_dispersion.mp4';
CFG.transverse_output_name = 'directional_cpp_large_spreading_eta11_transverse_linear_peak.mp4';
CFG.overwrite_existing = false;
CFG.phi_shifts_deg = [0 90 180 270];
CFG.lambda = 225;
CFG.gravity = 9.81;
CFG.max_frames = inf;
CFG.frame_stride = 1;
CFG.frame_rate = 8;
CFG.use_black_white_contours = true;
CFG.use_filled_contours = false;
CFG.contour_line_color = [0.35 0.35 0.35];
CFG.contour_line_width = 0.7;
CFG.num_contours_main = 8;
CFG.num_contours_subharmonic = 9;
CFG.center_coordinates = true;
CFG.normalize_by_A = false;
CFG.Akp = 0.15;
CFG.kp_ref = 2 * pi / CFG.lambda;
CFG.subharmonic_cutoff_factor = 1.20;
CFG.subharmonic_transition_factor = 0.35;
CFG.show_colorbar = false;
CFG.figure_position = [80 60 1380 980];
CFG.make_centerline_video = true;
CFG.make_transverse_video = true;
CFG.centerline_figure_position = [120 80 1380 900];
CFG.centerline_line_width = 1.6;
CFG.centerline_actual_color = [0.05 0.25 0.65];
CFG.centerline_linear_color = [0.80 0.20 0.10];
CFG.centerline_error_color = [0.20 0.20 0.20];
CFG.centerline_show_legend = true;

if ~isfolder(CFG.data_root)
    error('Data root not found: %s', CFG.data_root);
end

case_folders = resolve_phase_folders_local(CFG.data_root, CFG.phi_shifts_deg);
common_steps = intersect_ep_steps_local(case_folders);
selected_steps = common_steps(1:CFG.frame_stride:end);
if isfinite(CFG.max_frames)
    selected_steps = selected_steps(1:min(numel(selected_steps), CFG.max_frames));
end

if isempty(selected_steps)
    error('No common EP_*.bin snapshots found across the four phase folders.');
end

fprintf('Using %d common EP frames from %s\n', numel(selected_steps), CFG.data_root);

four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

sample_path = fullfile(case_folders{1}, sprintf('EP_%05d.bin', selected_steps(1)));
[X, Y, ~, ~] = read_ow3d_snapshot_local(sample_path);
x_vec = X(:, 1);
y_vec = Y(1, :);
kp = 2 * pi / CFG.lambda;
subharmonic_cutoff = CFG.subharmonic_cutoff_factor * kp;
subharmonic_transition = CFG.subharmonic_transition_factor * kp;
time_meta = read_focus_time_metadata_local(case_folders{1});
centerline_idx = resolve_centerline_index_local(y_vec);

if CFG.center_coordinates
    x_plot = x_vec - 0.5 * (x_vec(1) + x_vec(end));
    y_plot = y_vec - 0.5 * (y_vec(1) + y_vec(end));
else
    x_plot = x_vec;
    y_plot = y_vec;
end

scale_factor = 1.0;
scale_label = 'm';
if CFG.normalize_by_A
    A = CFG.Akp / CFG.kp_ref;
    if ~(isfinite(A) && A > 0)
        error('Invalid normalization scale derived from Akp and kp_ref.');
    end
    scale_factor = A;
    scale_label = 'scaled';
end

if ~isfolder(CFG.output_dir)
    mkdir(CFG.output_dir);
end

video_path = resolve_output_video_path_local(CFG.output_dir, CFG.output_name, CFG.overwrite_existing);
writer = VideoWriter(video_path, 'MPEG-4');
writer.FrameRate = CFG.frame_rate;
open(writer);

fig = figure('Color', 'w', 'Position', CFG.figure_position);
tile = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

if CFG.make_centerline_video || CFG.make_transverse_video
    t_rel_tp_series = zeros(numel(selected_steps), 1);
    t_case_seconds_series = zeros(numel(selected_steps), 1);
    linear_state = [];
    centerline_reference_time_seconds = NaN;
end

if CFG.make_centerline_video
    eta11_centerline_ow3d = zeros(numel(selected_steps), numel(x_vec));
    eta11_centerline_linear = zeros(numel(selected_steps), numel(x_vec));
end

if CFG.make_transverse_video
    eta11_transverse_ow3d = zeros(numel(selected_steps), numel(y_vec));
    eta11_transverse_linear = zeros(numel(selected_steps), numel(y_vec));
    transverse_x_peak_plot = zeros(numel(selected_steps), 1);
end

for frame_idx = 1:numel(selected_steps)
    ep_step = selected_steps(frame_idx);
    eta_phases = zeros(numel(CFG.phi_shifts_deg), numel(x_vec), numel(y_vec));

    for phase_idx = 1:numel(CFG.phi_shifts_deg)
        bin_path = fullfile(case_folders{phase_idx}, sprintf('EP_%05d.bin', ep_step));
        [X_cur, Y_cur, eta_cur, ~] = read_ow3d_snapshot_local(bin_path);
        if frame_idx == 1
            X = X_cur;
            Y = Y_cur;
        end
        eta_phases(phase_idx, :, :) = eta_cur;
    end

    eta_harmonics = reconstruct_harmonics_local(eta_phases, four_phase_coef);

    eta11 = frequency_filtering_2d_local(squeeze(eta_harmonics(1, :, :)), x_vec, y_vec, kp, 1);
    eta22 = frequency_filtering_2d_local(squeeze(eta_harmonics(2, :, :)), x_vec, y_vec, kp, 2);
    eta33 = frequency_filtering_2d_local(squeeze(eta_harmonics(3, :, :)), x_vec, y_vec, kp, 3);
    eta20 = lowpass_component_2d_local(squeeze(eta_harmonics(4, :, :)), x_vec, y_vec, subharmonic_cutoff, subharmonic_transition);
    t_rel_tp = compute_relative_focus_time_tp_local(ep_step, time_meta);
    t_case_seconds = compute_case_time_seconds_local(ep_step, time_meta);

    if CFG.make_centerline_video || CFG.make_transverse_video
        if isempty(linear_state)
            linear_state = build_linear_dispersion_state_2d_local(eta11, x_vec, y_vec, time_meta.depth, CFG.gravity);
            centerline_reference_time_seconds = t_case_seconds;
        end

        eta11_linear = propagate_linear_dispersion_field_2d_local(linear_state, t_case_seconds - centerline_reference_time_seconds);
        t_rel_tp_series(frame_idx) = t_rel_tp;
        t_case_seconds_series(frame_idx) = t_case_seconds;
    end

    if CFG.make_centerline_video
        eta11_centerline_ow3d(frame_idx, :) = eta11(:, centerline_idx).';
        eta11_centerline_linear(frame_idx, :) = eta11_linear(:, centerline_idx).';
    end

    if CFG.make_transverse_video
        linear_centerline_envelope = abs(hilbert(eta11_linear(:, centerline_idx)));
        [~, peak_x_idx] = max(linear_centerline_envelope);
        eta11_transverse_ow3d(frame_idx, :) = eta11(peak_x_idx, :);
        eta11_transverse_linear(frame_idx, :) = eta11_linear(peak_x_idx, :);
        transverse_x_peak_plot(frame_idx) = x_plot(peak_x_idx);
    end

    fields = {
        eta11 / scale_factor, ...
        eta22 / scale_factor, ...
        eta33 / scale_factor, ...
        eta20 / scale_factor};
    titles = {
        'First Order', ...
        'Second Order', ...
        'Third Order', ...
        'Second-Order Subharmonic'};

    clf(fig);
    tile = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('Directional OW3D harmonic evolution | EP %05d | t/T_p = %.2f relative to linear focus | frame %d/%d', ...
        ep_step, t_rel_tp, frame_idx, numel(selected_steps)), 'FontSize', 16, 'FontWeight', 'bold');

    for panel_idx = 1:4
        ax = nexttile(tile);
        draw_field_panel_local(ax, x_plot, y_plot, fields{panel_idx}, titles{panel_idx}, CFG, scale_label, panel_idx);
    end

    drawnow;
    writeVideo(writer, getframe(fig));

    if mod(frame_idx, 10) == 0 || frame_idx == numel(selected_steps)
        fprintf('Wrote frame %d / %d (EP_%05d)\n', frame_idx, numel(selected_steps), ep_step);
    end
end

close(writer);
fprintf('Saved video to %s\n', video_path);

if CFG.make_centerline_video
    centerline_video_path = resolve_output_video_path_local(CFG.output_dir, CFG.centerline_output_name, CFG.overwrite_existing);
    write_eta11_centerline_video_local(centerline_video_path, x_plot, y_plot(centerline_idx), selected_steps, ...
        t_rel_tp_series, eta11_centerline_ow3d / scale_factor, eta11_centerline_linear / scale_factor, CFG, scale_label);
    fprintf('Saved centerline eta11 video to %s\n', centerline_video_path);
end

if CFG.make_transverse_video
    transverse_video_path = resolve_output_video_path_local(CFG.output_dir, CFG.transverse_output_name, CFG.overwrite_existing);
    write_eta11_transverse_video_local(transverse_video_path, y_plot, selected_steps, t_rel_tp_series, ...
        transverse_x_peak_plot, eta11_transverse_ow3d / scale_factor, eta11_transverse_linear / scale_factor, CFG, scale_label);
    fprintf('Saved transverse eta11 video to %s\n', transverse_video_path);
end

function video_path = resolve_output_video_path_local(output_dir, output_name, overwrite_existing)
    video_path = fullfile(output_dir, output_name);
    if overwrite_existing
        if isfile(video_path)
            delete(video_path);
        end
        return;
    end

    if ~isfile(video_path)
        return;
    end

    [~, stem, ext] = fileparts(output_name);
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    video_path = fullfile(output_dir, sprintf('%s_%s%s', stem, timestamp, ext));
end

function meta = read_focus_time_metadata_local(case_folder)
    readme_path = fullfile(case_folder, 'OW_readme.txt');
    log_path = fullfile(case_folder, 'LOG.txt');

    if ~isfile(readme_path)
        error('Cannot find OW_readme.txt in %s', case_folder);
    end
    if ~isfile(log_path)
        error('Cannot find LOG.txt in %s', case_folder);
    end

    readme_text = fileread(readme_path);
    log_text = fileread(log_path);

    tp_tok = regexp(readme_text, 'Tp=(?<Tp>[-+]?\d*\.?\d+)', 'names', 'once');
    tinit_tok = regexp(readme_text, 'Initial condition time relative to focus:\s*(?<tinit>[-+]?\d*\.?\d+)\s*Tp', 'names', 'once');
    tend_tok = regexp(readme_text, 'Target end time relative to focus:\s*(?<tend>[-+]?\d*\.?\d+)\s*Tp', 'names', 'once');
    stride_tok = regexp(readme_text, 'Surface output stride:\s*every\s*(?<stride>\d+)\s*time step', 'names', 'once');
    nsteps_tok = regexp(log_text, 'Number of time steps chosen:\s*(?<nsteps>\d+)', 'names', 'once');
    dt_tok = regexp(log_text, 'Size of time increment:\s*(?<dt>[.\deE+-]+)', 'names', 'once');
    depth_tok = regexp(readme_text, 'h=(?<depth>[-+]?\d*\.?\d+)', 'names', 'once');

    if isempty(tp_tok) || isempty(tinit_tok) || isempty(tend_tok) || isempty(stride_tok) || ...
            isempty(nsteps_tok) || isempty(dt_tok) || isempty(depth_tok)
        error('Could not parse full timing/depth metadata from %s and %s', readme_path, log_path);
    end

    meta = struct();
    meta.Tp = str2double(tp_tok.Tp);
    meta.t_init_rel_focus_tp = str2double(tinit_tok.tinit);
    meta.t_end_rel_focus_tp = str2double(tend_tok.tend);
    meta.surface_stride = str2double(stride_tok.stride);
    meta.n_steps = str2double(nsteps_tok.nsteps);
    meta.dt = str2double(dt_tok.dt);
    meta.depth = str2double(depth_tok.depth);
end

function t_rel_tp = compute_relative_focus_time_tp_local(ep_step, meta)
    if ep_step == 99999
        t_rel_tp = meta.t_end_rel_focus_tp;
        return;
    end

    elapsed_time = ep_step * meta.dt;
    t_rel_tp = meta.t_init_rel_focus_tp + elapsed_time / meta.Tp;
end

function t_case_seconds = compute_case_time_seconds_local(ep_step, meta)
    if ep_step == 99999
        t_case_seconds = (meta.t_end_rel_focus_tp - meta.t_init_rel_focus_tp) * meta.Tp;
        return;
    end

    t_case_seconds = ep_step * meta.dt;
end

function case_folders = resolve_phase_folders_local(data_root, phi_shifts_deg)
    listing = dir(fullfile(data_root, '*phi_*'));
    listing = listing([listing.isdir]);

    case_folders = cell(1, numel(phi_shifts_deg));
    for i = 1:numel(phi_shifts_deg)
        pattern = sprintf('phi_%d', phi_shifts_deg(i));
        match_idx = find(contains({listing.name}, pattern), 1, 'first');
        if isempty(match_idx)
            error('Could not find folder for phase shift %d deg under %s', phi_shifts_deg(i), data_root);
        end
        case_folders{i} = fullfile(data_root, listing(match_idx).name);
    end
end

function common_steps = intersect_ep_steps_local(case_folders)
    step_sets = cell(1, numel(case_folders));
    for i = 1:numel(case_folders)
        files = dir(fullfile(case_folders{i}, 'EP_*.bin'));
        steps = nan(1, numel(files));
        for j = 1:numel(files)
            token = regexp(files(j).name, '^EP_(\d+)\.bin$', 'tokens', 'once');
            if ~isempty(token)
                steps(j) = str2double(token{1});
            end
        end
        step_sets{i} = steps(isfinite(steps));
    end

    common_steps = step_sets{1};
    for i = 2:numel(step_sets)
        common_steps = intersect(common_steps, step_sets{i});
    end
    common_steps = sort(common_steps(:).');
end

function [X, Y, eta_field, phi_field] = read_ow3d_snapshot_local(bin_path)
    byteorder = 'ieee-le';
    fid = fopen(bin_path, 'r', byteorder);
    if fid < 0
        error('Unable to open OW3D bin file: %s', bin_path);
    end
    cleanup = onCleanup(@() fclose(fid));

    fread(fid, 1, 'int32');
    nx = fread(fid, 1, 'int32');
    ny = fread(fid, 1, 'int32');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    X = fread(fid, [nx ny], 'float64');
    Y = fread(fid, [nx ny], 'float64');
    fread(fid, 1, 'int32');

    fread(fid, 1, 'int32');
    eta_field = fread(fid, [nx ny], 'float64');
    phi_field = fread(fid, [nx ny], 'float64');
end

function eta_harmonics = reconstruct_harmonics_local(eta_phases, coef)
    eta_hilbert = hilbert2d_all_local(eta_phases);
    all_eta = cat(1, real(eta_phases), -imag(eta_hilbert));

    eta_harmonics = zeros(4, size(eta_phases, 2), size(eta_phases, 3));
    for n = 1:4
        weights = reshape(coef(n, :), [8 1 1]);
        eta_harmonics(n, :, :) = sum(all_eta .* weights, 1);
    end
end

function X_hilbert = hilbert2d_all_local(X)
    [n_phase, nx, ny] = size(X);
    X_hilbert = zeros(size(X));

    for p = 1:n_phase
        FX = fft2(squeeze(X(p, :, :)));
        H = zeros(nx, ny);
        H(1:floor(nx / 2), :) = 1;
        H(ceil(nx / 2):end, :) = -1;
        X_hilbert(p, :, :) = ifft2(FX .* H);
    end
end

function field_out = frequency_filtering_2d_local(field_in, x_vec, y_vec, kp, n)
    x_vec = x_vec(:);
    y_vec = y_vec(:);
    [nx, ny] = size(field_in);

    dx = x_vec(2) - x_vec(1);
    dy = y_vec(2) - y_vec(1);
    dkx = 2 * pi / (nx * dx);
    dky = 2 * pi / (ny * dy);

    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    ky = [0:ceil(ny / 2) - 1, -floor(ny / 2):-1] * dky;
    [KX, KY] = ndgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);

    sigma = 0.45 * kp;
    k_target = n * kp;
    mask = exp(-((K - k_target).^2) / (2 * sigma^2));

    field_fft = fft2(field_in);
    field_out = real(ifft2(field_fft .* mask));
end

function field_out = lowpass_component_2d_local(field_in, x_vec, y_vec, k_cutoff, transition)
    x_vec = x_vec(:);
    y_vec = y_vec(:);
    [nx, ny] = size(field_in);

    dx = x_vec(2) - x_vec(1);
    dy = y_vec(2) - y_vec(1);
    dkx = 2 * pi / (nx * dx);
    dky = 2 * pi / (ny * dy);

    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    ky = [0:ceil(ny / 2) - 1, -floor(ny / 2):-1] * dky;
    [KX, KY] = ndgrid(kx, ky);
    K = sqrt(KX.^2 + KY.^2);

    mask = exp(-(K / max(transition, min(dkx, dky))).^4);
    mask(K <= k_cutoff) = 1;

    field_fft = fft2(field_in);
    field_out = real(ifft2(field_fft .* mask));
end

function idx = resolve_centerline_index_local(y_vec)
    y_vec = y_vec(:);
    [~, idx] = min(abs(y_vec - 0.5 * (y_vec(1) + y_vec(end))));
end

function state = build_linear_dispersion_state_2d_local(field0, x_vec, y_vec, depth, gravity)
    x_vec = x_vec(:);
    y_vec = y_vec(:);
    [nx, ny] = size(field0);

    dx = x_vec(2) - x_vec(1);
    dy = y_vec(2) - y_vec(1);
    dkx = 2 * pi / (nx * dx);
    dky = 2 * pi / (ny * dy);

    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    ky = [0:ceil(ny / 2) - 1, -floor(ny / 2):-1] * dky;
    [KX, KY] = ndgrid(kx, ky);
    K = hypot(KX, KY);

    omega = sqrt(gravity * K .* tanh(K * depth));
    omega(~isfinite(omega)) = 0;
    omega(K <= eps) = 0;

    state = struct();
    state.field_fft0 = fft2(field0);
    state.omega = omega;
    state.pos_mask = (KX > 0) | ((KX == 0) & (KY > 0));
    state.neg_mask = (KX < 0) | ((KX == 0) & (KY < 0));
end

function field_t = propagate_linear_dispersion_field_2d_local(state, delta_t)
    phase = ones(size(state.field_fft0));
    phase(state.pos_mask) = exp(-1i * state.omega(state.pos_mask) * delta_t);
    phase(state.neg_mask) = exp(1i * state.omega(state.neg_mask) * delta_t);
    field_t = real(ifft2(state.field_fft0 .* phase));
end

function write_eta11_centerline_video_local(video_path, x_plot, y_centerline, selected_steps, t_rel_tp_series, ...
        eta11_centerline_ow3d, eta11_centerline_linear, CFG, scale_label)
    eta_error = eta11_centerline_ow3d - eta11_centerline_linear;
    amp_peak = max(abs([eta11_centerline_ow3d(:); eta11_centerline_linear(:)]));
    err_peak = max(abs(eta_error(:)));

    if ~(isfinite(amp_peak) && amp_peak > 0)
        amp_peak = 1;
    end
    if ~(isfinite(err_peak) && err_peak > 0)
        err_peak = 1;
    end

    writer = VideoWriter(video_path, 'MPEG-4');
    writer.FrameRate = CFG.frame_rate;
    open(writer);

    fig = figure('Color', 'w', 'Position', CFG.centerline_figure_position);

    for frame_idx = 1:numel(selected_steps)
        clf(fig);
        tile = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tile, sprintf(['Centerline \\eta_{11} evolution | y = %.2f m | EP %05d | t/T_p = %.2f ' ...
            '| linear dispersion from first frame | frame %d/%d'], ...
            y_centerline, selected_steps(frame_idx), t_rel_tp_series(frame_idx), frame_idx, numel(selected_steps)), ...
            'FontSize', 16, 'FontWeight', 'bold');

        ax1 = nexttile(tile);
        plot(ax1, x_plot, eta11_centerline_ow3d(frame_idx, :), 'Color', CFG.centerline_actual_color, ...
            'LineWidth', CFG.centerline_line_width);
        hold(ax1, 'on');
        plot(ax1, x_plot, eta11_centerline_linear(frame_idx, :), '--', 'Color', CFG.centerline_linear_color, ...
            'LineWidth', CFG.centerline_line_width);
        grid(ax1, 'on');
        box(ax1, 'on');
        ylabel(ax1, sprintf('\\eta_{11} (%s)', scale_label));
        xlim(ax1, [x_plot(1), x_plot(end)]);
        ylim(ax1, 1.08 * amp_peak * [-1 1]);
        title(ax1, 'OW3D vs linearly dispersed centerline', 'FontSize', 13, 'FontWeight', 'bold');
        if CFG.centerline_show_legend
            legend(ax1, {'OW3D centerline', 'Linear dispersion'}, 'Location', 'northwest');
        end

        ax2 = nexttile(tile);
        plot(ax2, x_plot, eta_error(frame_idx, :), 'Color', CFG.centerline_error_color, ...
            'LineWidth', CFG.centerline_line_width);
        grid(ax2, 'on');
        box(ax2, 'on');
        xlabel(ax2, 'x (m)');
        ylabel(ax2, sprintf('\\Delta\\eta_{11} (%s)', scale_label));
        xlim(ax2, [x_plot(1), x_plot(end)]);
        ylim(ax2, 1.08 * err_peak * [-1 1]);
        title(ax2, '\eta_{11}^{OW3D} - \eta_{11}^{linear}', 'FontSize', 13, 'FontWeight', 'bold');

        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function write_eta11_transverse_video_local(video_path, y_plot, selected_steps, t_rel_tp_series, x_peak_plot, ...
        eta11_transverse_ow3d, eta11_transverse_linear, CFG, scale_label)
    eta_error = eta11_transverse_ow3d - eta11_transverse_linear;
    amp_peak = max(abs([eta11_transverse_ow3d(:); eta11_transverse_linear(:)]));
    err_peak = max(abs(eta_error(:)));

    if ~(isfinite(amp_peak) && amp_peak > 0)
        amp_peak = 1;
    end
    if ~(isfinite(err_peak) && err_peak > 0)
        err_peak = 1;
    end

    writer = VideoWriter(video_path, 'MPEG-4');
    writer.FrameRate = CFG.frame_rate;
    open(writer);

    fig = figure('Color', 'w', 'Position', CFG.centerline_figure_position);

    for frame_idx = 1:numel(selected_steps)
        clf(fig);
        tile = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
        title(tile, sprintf(['Transverse \\eta_{11} evolution | x = %.2f m (linear-envelope peak) | EP %05d ' ...
            '| t/T_p = %.2f | frame %d/%d'], ...
            x_peak_plot(frame_idx), selected_steps(frame_idx), t_rel_tp_series(frame_idx), frame_idx, numel(selected_steps)), ...
            'FontSize', 16, 'FontWeight', 'bold');

        ax1 = nexttile(tile);
        plot(ax1, y_plot, eta11_transverse_ow3d(frame_idx, :), 'Color', CFG.centerline_actual_color, ...
            'LineWidth', CFG.centerline_line_width);
        hold(ax1, 'on');
        plot(ax1, y_plot, eta11_transverse_linear(frame_idx, :), '--', 'Color', CFG.centerline_linear_color, ...
            'LineWidth', CFG.centerline_line_width);
        grid(ax1, 'on');
        box(ax1, 'on');
        ylabel(ax1, sprintf('\\eta_{11} (%s)', scale_label));
        xlim(ax1, [y_plot(1), y_plot(end)]);
        ylim(ax1, 1.08 * amp_peak * [-1 1]);
        title(ax1, 'OW3D vs linearly dispersed transverse slice', 'FontSize', 13, 'FontWeight', 'bold');
        if CFG.centerline_show_legend
            legend(ax1, {'OW3D transverse slice', 'Linear dispersion'}, 'Location', 'northwest');
        end

        ax2 = nexttile(tile);
        plot(ax2, y_plot, eta_error(frame_idx, :), 'Color', CFG.centerline_error_color, ...
            'LineWidth', CFG.centerline_line_width);
        grid(ax2, 'on');
        box(ax2, 'on');
        xlabel(ax2, 'y (m)');
        ylabel(ax2, sprintf('\\Delta\\eta_{11} (%s)', scale_label));
        xlim(ax2, [y_plot(1), y_plot(end)]);
        ylim(ax2, 1.08 * err_peak * [-1 1]);
        title(ax2, '\eta_{11}^{OW3D} - \eta_{11}^{linear}', 'FontSize', 13, 'FontWeight', 'bold');

        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function draw_field_panel_local(ax, x_plot, y_plot, field_in, title_text, CFG, scale_label, panel_idx)
    field_in = double(field_in);
    peak_abs = max(abs(field_in(:)));
    field_span = max(field_in(:)) - min(field_in(:));

    if panel_idx == 4
        num_contours = CFG.num_contours_subharmonic;
    else
        num_contours = CFG.num_contours_main;
    end

    if peak_abs < eps || field_span < 1e-12
        imagesc(ax, x_plot, y_plot, field_in.');
        colormap(ax, gray(256));
        caxis(ax, [-1 1]);
    elseif CFG.use_black_white_contours
        levels = linspace(-peak_abs, peak_abs, num_contours);
        if CFG.use_filled_contours
            contourf(ax, x_plot, y_plot, field_in.', levels, 'LineColor', CFG.contour_line_color, ...
                'LineWidth', CFG.contour_line_width);
            colormap(ax, gray(256));
            caxis(ax, [-peak_abs peak_abs]);
        else
            contour(ax, x_plot, y_plot, field_in.', levels, 'LineColor', CFG.contour_line_color, ...
                'LineWidth', CFG.contour_line_width);
        end
    else
        imagesc(ax, x_plot, y_plot, field_in.');
        colormap(ax, parula(256));
    end

    axis(ax, 'image');
    set(ax, 'YDir', 'normal', 'FontSize', 12, 'LineWidth', 1.0);
    xlabel(ax, 'x (m)');
    ylabel(ax, 'y (m)');
    title(ax, sprintf('%s | max = %.3e %s', title_text, max(abs(field_in(:))), scale_label), ...
        'FontSize', 13, 'FontWeight', 'bold');
    box(ax, 'on');

    if CFG.show_colorbar
        colorbar(ax);
    end
end
