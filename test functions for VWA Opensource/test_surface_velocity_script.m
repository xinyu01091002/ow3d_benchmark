%% Test Script for VWA Surface Velocity Functions (u and w)
% Compares OW3D surface harmonic reconstructions against the reusable VWA
% helpers driven only by the linear surface elevation eta^(1).

clc; clear; close all;

%% 1. Paths and configuration
current_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(current_dir);
raw_dir = fullfile(project_dir, 'processed_boundkinematics');
export_folder = fullfile(current_dir, 'figures_comparison');
if ~isfolder(export_folder), mkdir(export_folder); end

config = struct();
config.snapshot_tag = '';
config.lambda = 225;
config.gravity = 9.81;
config.small_kd_min = 0.3;
config.apply_eta11_filter = false;
config.sigma_mode = 'surface'; % 'surface' or 'index'
config.sigma_index = [];
config.y_index = 1;

snapshot_path = resolve_snapshot_path(raw_dir, config.snapshot_tag);
fprintf('Loading raw kinematics snapshot from: %s\n', snapshot_path);
S = load(snapshot_path);

raw_snapshot = S.raw_snapshot;
raw_meta = S.raw_meta;
phase_names = fieldnames(raw_snapshot);
phase_names = sort(phase_names);
if numel(phase_names) ~= 4
    error('Expected four stored phases, found %d.', numel(phase_names));
end

x_vec = raw_meta.x(:);
sigma_vec = raw_meta.sigma(:);
depth_value = mean(raw_snapshot.(phase_names{1}).h(:), 'omitnan');
kp = 2 * pi / config.lambda;

if strcmpi(config.sigma_mode, 'surface')
    [~, sigma_idx] = max(sigma_vec);
else
    sigma_idx = config.sigma_index;
end
if isempty(sigma_idx) || sigma_idx < 1 || sigma_idx > numel(sigma_vec)
    error('Invalid sigma index for the surface-velocity test script.');
end

eta_phases = zeros(numel(phase_names), numel(x_vec));
u_phases = zeros(numel(phase_names), numel(x_vec));
w_phases = zeros(numel(phase_names), numel(x_vec));

for idx = 1:numel(phase_names)
    snapshot = raw_snapshot.(phase_names{idx});
    eta_phases(idx, :) = squeeze(snapshot.eta(:, config.y_index));
    u_slice = squeeze(snapshot.u(sigma_idx, :, config.y_index));
    w_slice = squeeze(snapshot.w(sigma_idx, :, config.y_index));
    u_phases(idx, :) = u_slice(:).';
    w_phases(idx, :) = w_slice(:).';
end

coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

eta_h = reconstruct_harmonics_1d_local(eta_phases, coef);
u_h = reconstruct_harmonics_1d_local(u_phases, coef);
w_h = reconstruct_harmonics_1d_local(w_phases, coef);

if config.apply_eta11_filter
    eta11 = frequency_filtering_1d_local(eta_h(1, :), x_vec, kp, 1).';
else
    eta11 = eta_h(1, :).';
end

opts = struct('analytic_side', 'neg', 'small_kd_min', config.small_kd_min);
u_vwa = vwa_compute_surface_quantity(eta11, x_vec, depth_value, config.gravity, 'u', opts);
w_vwa = vwa_compute_surface_quantity(eta11, x_vec, depth_value, config.gravity, 'w', opts);

x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / config.lambda;
export_surface_velocity_compare(export_folder, x_plot, 'u', u_h, u_vwa, sigma_vec(sigma_idx), raw_meta.time_value);
export_surface_velocity_compare(export_folder, x_plot, 'w', w_h, w_vwa, sigma_vec(sigma_idx), raw_meta.time_value);

disp('Surface velocity VWA comparison complete.');

%% Local helper functions
function h = reconstruct_harmonics_1d_local(phase_data, coef)
    n_phase = size(phase_data, 1);
    nx = size(phase_data, 2);
    if n_phase ~= 4
        error('Need exactly four phase-shifted inputs for the harmonic reconstruction.');
    end

    h = zeros(4, nx);
    for ix = 1:nx
        y = phase_data(:, ix);
        analytic_part = imag(hilbert(y));
        state = [y(1); analytic_part(1); y(2); analytic_part(2); y(3); analytic_part(3); y(4); analytic_part(4)];
        h(:, ix) = coef * state;
    end
end

function eta_filtered = frequency_filtering_1d_local(eta, x, kp, order)
    eta = eta(:).';
    x = x(:).';
    nx = numel(x);
    dx = mean(diff(x));
    dk = 2 * pi / (nx * dx);
    k = fftshift((-floor(nx / 2):ceil(nx / 2) - 1) * dk);
    band = 0.35 * kp;
    spectrum = fftshift(fft(eta));
    mask = abs(abs(k) - order * kp) <= band;
    spectrum(~mask) = 0;
    eta_filtered = real(ifft(ifftshift(spectrum)));
end

function export_surface_velocity_compare(export_folder, x_plot, quantity_name, ow3d_h, vwa_out, sigma_value, time_value)
    fig = figure('Color', 'w', 'Position', [80 80 1450 920]);
    tile = tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('Surface %s: OW3D vs VWA (\\sigma = %.3f, t = %.4f s)', ...
        quantity_name, sigma_value, time_value), ...
        'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    order_titles = {'(a) First harmonic', '(b) Second harmonic', '(c) Third harmonic'};
    vwa_h = {vwa_out.order1(:), vwa_out.order2(:), vwa_out.order3(:)};
    diagnostics = strings(3, 1);

    for order = 1:3
        ow3d_vals = ow3d_h(order, :).';
        vwa_vals = vwa_h{order};
        y_limits = paired_ylim_local(ow3d_vals, vwa_vals);

        ax = nexttile(tile);
        plot(ax, x_plot, ow3d_vals, 'k-', 'LineWidth', 1.6, 'DisplayName', 'OW3D'); hold(ax, 'on');
        plot(ax, x_plot, vwa_vals, '--', 'Color', [0.80 0.26 0.18], 'LineWidth', 1.6, 'DisplayName', 'VWA');
        grid(ax, 'on'); box(ax, 'on');
        ylabel(ax, quantity_axis_label_local(quantity_name), 'Interpreter', 'latex', 'FontSize', 12);
        title(ax, order_titles{order}, 'Interpreter', 'tex', 'FontSize', 13);
        set(ax, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
        xlim(ax, [min(x_plot), max(x_plot)]);
        ylim(ax, y_limits);
        if order == 1
            legend(ax, 'Location', 'best');
        end
        if order == 3
            xlabel(ax, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 12);
        end

        metrics = compare_metrics_local(ow3d_vals, vwa_vals);
        diagnostics(order) = sprintf('order %d: corr = %.3f, RMSE = %.3e, peak ratio = %.3f', ...
            order, metrics.corr, metrics.rmse, metrics.peak_ratio);
        text(ax, 0.02, 0.92, diagnostics(order), 'Units', 'normalized', ...
            'VerticalAlignment', 'top', 'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);
    end

    annotation(fig, 'textbox', [0.13 0.01 0.82 0.05], ...
        'String', strjoin(cellstr(diagnostics), '    |    '), ...
        'Interpreter', 'none', 'EdgeColor', 'none', 'HorizontalAlignment', 'left', ...
        'FontName', 'Times New Roman', 'FontSize', 10);

    exportgraphics(fig, fullfile(export_folder, sprintf('comparison_surface_%s.png', lower(quantity_name))), ...
        'Resolution', 300);
end

function metrics = compare_metrics_local(reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    metrics = struct();
    if all(abs(reference) < eps) || all(abs(candidate) < eps)
        metrics.corr = NaN;
    else
        cc = corrcoef(reference, candidate);
        metrics.corr = cc(1, 2);
    end
    metrics.rmse = sqrt(mean((reference - candidate).^2));
    ref_peak = max(abs(reference));
    cand_peak = max(abs(candidate));
    if ref_peak <= eps
        metrics.peak_ratio = NaN;
    else
        metrics.peak_ratio = cand_peak / ref_peak;
    end
end

function y_limits = paired_ylim_local(a, b)
    max_abs = max(abs([a(:); b(:)]));
    max_abs = max(max_abs, 1e-8);
    y_limits = [-1, 1] * 1.10 * max_abs;
end

function label = quantity_axis_label_local(quantity_name)
    switch lower(quantity_name)
        case 'u'
            label = '$u_s$ (m/s)';
        case 'w'
            label = '$w_s$ (m/s)';
        otherwise
            label = ['$', quantity_name, '$'];
    end
end

function snapshot_path = resolve_snapshot_path(raw_dir, snapshot_tag)
    if ~isempty(snapshot_tag)
        snapshot_path = fullfile(raw_dir, snapshot_tag);
        if ~isfile(snapshot_path)
            error('Configured raw snapshot not found: %s', snapshot_path);
        end
        return;
    end

    candidates = dir(fullfile(raw_dir, 'OW3D_boundkinematics_raw_tidx_*.mat'));
    if isempty(candidates)
        error('No raw bound-kinematics snapshots found under %s.', raw_dir);
    end

    [~, idx] = max([candidates.datenum]);
    snapshot_path = fullfile(candidates(idx).folder, candidates(idx).name);
end
