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
config.compare_mf12_subharmonic = true;
config.subharmonic_cutoff_factor = 1.0;
config.subharmonic_transition_factor = 0.35;
config.mf12_linear_energy_keep = 0.99;

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
    eta_slice = squeeze(snapshot.eta(config.y_index, :));
    u_slice = squeeze(snapshot.u(sigma_idx, :));
    w_slice = squeeze(snapshot.w(sigma_idx, :));
    eta_phases(idx, :) = eta_slice(:).';
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

if config.compare_mf12_subharmonic
    addpath(fullfile(project_dir, 'irregularWavesMF12', 'Source'));

    subharmonic_cutoff = config.subharmonic_cutoff_factor * kp;
    ow3d_u20 = lowpass_wavenumber_component_local(mean(u_phases, 1).', x_vec, subharmonic_cutoff, ...
        config.subharmonic_transition_factor * subharmonic_cutoff);
    ow3d_w20 = lowpass_wavenumber_component_local(mean(w_phases, 1).', x_vec, subharmonic_cutoff, ...
        config.subharmonic_transition_factor * subharmonic_cutoff);

    eta_linear_for_mf12 = eta_h(1, :).';
    mf12_sub = compute_mf12_second_order_filtered_surface(eta_linear_for_mf12, x_vec, depth_value, ...
        config.gravity, subharmonic_cutoff, config.subharmonic_transition_factor * subharmonic_cutoff, ...
        config.mf12_linear_energy_keep);

    export_subharmonic_compare(export_folder, x_plot, 'u', ow3d_u20, mf12_sub.u20, sigma_vec(sigma_idx), ...
        raw_meta.time_value, subharmonic_cutoff / kp, 'MF12');
    export_subharmonic_compare(export_folder, x_plot, 'w', ow3d_w20, mf12_sub.w20, sigma_vec(sigma_idx), ...
        raw_meta.time_value, subharmonic_cutoff / kp, 'MF12');
end

disp('Surface velocity VWA comparison complete.');

%% Local helper functions
function h = reconstruct_harmonics_1d_local(phase_data, coef)
    n_phase = size(phase_data, 1);
    nx = size(phase_data, 2);
    if n_phase ~= 4
        error('Need exactly four phase-shifted inputs for the harmonic reconstruction.');
    end

    analytic_part = hilbert(phase_data.').';
    all_fields = zeros(8, nx);
    all_fields(1:2:end, :) = real(phase_data);
    all_fields(2:2:end, :) = -imag(analytic_part);
    h = zeros(4, nx);
    for n = 1:4
        h(n, :) = coef(n, :) * all_fields;
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

function out = compute_mf12_second_order_filtered_surface(eta11, x_vec, h, g, k_cutoff, transition, energy_keep)
    eta11 = eta11(:);
    x_vec = x_vec(:);
    nx = numel(x_vec);
    dx = x_vec(2) - x_vec(1);
    kx_grid = vwa_kxgrid(nx, dx);
    fft_eta = fft(eta11) / nx;

    if mod(nx, 2) == 0
        positive_idx = 2:(nx / 2);
    else
        positive_idx = 2:((nx + 1) / 2);
    end

    positive_idx = positive_idx(abs(fft_eta(positive_idx)) > 1e-12 * max(abs(fft_eta)));
    positive_idx = select_energy_dominant_indices_local(fft_eta, positive_idx, energy_keep);
    kx = kx_grid(positive_idx).';
    ky = zeros(size(kx));
    a = 2 * real(fft_eta(positive_idx)).';
    b = 2 * imag(fft_eta(positive_idx)).';

    coeffs1 = mf12_direct_coefficients(1, g, h, a, b, kx, ky, 0, 0, 0);
    coeffs2 = mf12_direct_coefficients(2, g, h, a, b, kx, ky, 0, 0, 0);

    [u1, ~, w1, ~, phi1] = kinematicsMF12(1, coeffs1, x_vec.', 0, 0, 0); %#ok<ASGLU>
    [u2, ~, w2, ~, phi2] = kinematicsMF12(2, coeffs2, x_vec.', 0, 0, 0); %#ok<ASGLU>

    theta = coeffs1.omega(:) .* 0 - coeffs1.kx(:) .* x_vec.' - coeffs1.ky(:) .* 0;
    eta1_matrix = coeffs1.a(:) .* cos(theta) + coeffs1.b(:) .* sin(theta);
    eta1 = sum(eta1_matrix, 1).';

    theta2 = 2 * coeffs2.omega(:) .* 0 - coeffs2.kx_2(:) .* x_vec.' - coeffs2.ky_2(:) .* 0;
    eta2_self = sum(coeffs2.G_2(:) .* (coeffs2.A_2(:) .* cos(theta2) + coeffs2.B_2(:) .* sin(theta2)), 1).';

    eta2_pair = zeros(nx, 1);
    cnm = 0;
    for n = 1:coeffs2.N
        for m = (n + 1):coeffs2.N
            for pm = [1 -1]
                cnm = cnm + 1;
                theta_npm = coeffs2.omega_npm(cnm) .* 0 - coeffs2.kx_npm(cnm) .* x_vec - coeffs2.ky_npm(cnm) .* 0;
                eta2_pair = eta2_pair + coeffs2.G_npm(cnm) .* ...
                    (coeffs2.A_npm(cnm) .* cos(theta_npm) + coeffs2.B_npm(cnm) .* sin(theta_npm));
            end
        end
    end

    eta20 = (eta2_self + eta2_pair);
    phi20 = phi2(:) - phi1(:);
    u20 = u2(:) - u1(:);
    w20 = w2(:) - w1(:);

    eta20 = lowpass_wavenumber_component_local(eta20, x_vec, k_cutoff, transition);
    phi20 = lowpass_wavenumber_component_local(phi20, x_vec, k_cutoff, transition);
    u20 = lowpass_wavenumber_component_local(u20, x_vec, k_cutoff, transition);
    w20 = lowpass_wavenumber_component_local(w20, x_vec, k_cutoff, transition);

    out = struct();
    out.eta20 = eta20(:);
    out.phi20 = phi20(:);
    out.u20 = u20(:);
    out.w20 = w20(:);
    out.coeffs1 = coeffs1;
    out.coeffs2 = coeffs2;
    out.linear_indices = positive_idx(:);
    out.linear_energy_keep = energy_keep;
end

function field_out = lowpass_wavenumber_component_local(field_in, x_vec, k_cutoff, transition)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = exp(-(abs(kx) / max(transition, dkx)).^4);
    mask(abs(kx) <= k_cutoff) = 1;
    field_out = real(ifft(fft(field_in) .* mask));
end

function field_out = bandpass_harmonic_component_local(field_in, x_vec, k_center, k_sigma)
    field_in = field_in(:);
    x_vec = x_vec(:);
    nx = numel(field_in);
    dx = x_vec(2) - x_vec(1);
    dkx = 2 * pi / (nx * dx);
    kx = [0:ceil(nx / 2) - 1, -floor(nx / 2):-1]' * dkx;
    mask = exp(-((abs(kx) - k_center).^2) / (2 * k_sigma^2));
    field_out = real(ifft(fft(field_in) .* mask));
end

function keep_idx = select_energy_dominant_indices_local(fft_eta, candidate_idx, energy_keep)
    if nargin < 3 || isempty(energy_keep)
        energy_keep = 1.0;
    end
    energy_keep = min(max(energy_keep, 0), 1);
    if isempty(candidate_idx)
        keep_idx = candidate_idx;
        return;
    end

    spectral_energy = abs(fft_eta(candidate_idx)).^2;
    total_energy = sum(spectral_energy);
    if total_energy <= 0 || energy_keep >= 1
        keep_idx = candidate_idx;
        return;
    end

    [sorted_energy, order] = sort(spectral_energy, 'descend');
    cumulative_energy = cumsum(sorted_energy) / total_energy;
    cutoff_pos = find(cumulative_energy >= energy_keep, 1, 'first');
    keep_unsorted = candidate_idx(order(1:cutoff_pos));
    keep_idx = sort(keep_unsorted);
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

function export_subharmonic_compare(export_folder, x_plot, quantity_name, ow3d_sub, model_sub, sigma_value, time_value, cutoff_ratio, model_name)
    fig = figure('Color', 'w', 'Position', [80 80 1450 720]);
    tile = tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
    title(tile, sprintf('Second subharmonic surface %s: OW3D vs %s (\\sigma = %.3f, t = %.4f s, |k| < %.2f k_p)', ...
        quantity_name, model_name, sigma_value, time_value, cutoff_ratio), ...
        'Interpreter', 'tex', 'FontSize', 16, 'FontWeight', 'bold');

    ax1 = nexttile(tile);
    plot(ax1, x_plot, ow3d_sub, 'k-', 'LineWidth', 1.6, 'DisplayName', 'OW3D'); hold(ax1, 'on');
    plot(ax1, x_plot, model_sub, '--', 'Color', [0.80 0.26 0.18], 'LineWidth', 1.6, 'DisplayName', model_name);
    grid(ax1, 'on'); box(ax1, 'on');
    ylabel(ax1, quantity_axis_label_local(quantity_name), 'Interpreter', 'latex', 'FontSize', 12);
    title(ax1, 'Subharmonic comparison', 'Interpreter', 'tex', 'FontSize', 13);
    set(ax1, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlim(ax1, [min(x_plot), max(x_plot)]);
    ylim(ax1, paired_ylim_local(ow3d_sub, model_sub));
    legend(ax1, 'Location', 'best');

    metrics = compare_metrics_local(ow3d_sub, model_sub);
    text(ax1, 0.02, 0.92, sprintf('corr = %.3f, RMSE = %.3e, peak ratio = %.3f', ...
        metrics.corr, metrics.rmse, metrics.peak_ratio), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'FontSize', 10);

    ax2 = nexttile(tile);
    plot(ax2, x_plot, ow3d_sub - model_sub, '-', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.6);
    grid(ax2, 'on'); box(ax2, 'on');
    xlabel(ax2, '$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(ax2, quantity_axis_label_local(quantity_name), 'Interpreter', 'latex', 'FontSize', 12);
    title(ax2, sprintf('Difference (OW3D - %s)', model_name), 'Interpreter', 'tex', 'FontSize', 13);
    set(ax2, 'FontName', 'Times New Roman', 'FontSize', 12, 'LineWidth', 1.0);
    xlim(ax2, [min(x_plot), max(x_plot)]);
    ylim(ax2, paired_ylim_local(ow3d_sub - model_sub, zeros(size(ow3d_sub))));

    exportgraphics(fig, fullfile(export_folder, sprintf('comparison_%s_subharmonic_%s.png', lower(model_name), lower(quantity_name))), ...
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
