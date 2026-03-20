clc;
clear;
close all;

cfg = struct();
cfg.data_root = fullfile(pwd, 'uni initial condition', 'ow3d_kinematics_check');
cfg.folder_pattern = 'T_init-20_Tp_Alpha_1.0_Akp_006_kd1.0_phi_%d';
cfg.phi_shifts_deg = 0:90:270;
cfg.kinematics_file_id = 1;
cfg.time_index = 66;
cfg.lambda = 225;
cfg.gravity = 9.81;
cfg.apply_x_filter = true;

reader_dir = 'c:\Users\spet5947\OneDrive - Nexus365\Research\code transfer\function';
if isfolder(reader_dir)
    addpath(reader_dir);
else
    error('Reader directory not found: %s', reader_dir);
end

eta_phases = [];
phi_surface_phases = [];
phit_corr_surface_phases = [];
phit_uncorr_surface_phases = [];
sigma_vec = [];
x_vec = [];
t_selected = [];

for idx = 1:numel(cfg.phi_shifts_deg)
    case_folder = fullfile(cfg.data_root, sprintf(cfg.folder_pattern, cfg.phi_shifts_deg(idx)));
    kin_path = fullfile(case_folder, sprintf('Kinematics%02d.bin', cfg.kinematics_file_id));
    if ~isfile(kin_path)
        error('Missing kinematics file: %s', kin_path);
    end

    [it, eta, etat_m, ~, phi, phit_m, ~, ~, u, ~, w, uz, ~, ~, x, y, sigma, t] = Read_Kinamatics(kin_path, cfg.kinematics_file_id); %#ok<ASGLU>
    if cfg.time_index > it
        error('Requested time index %d exceeds available steps %d.', cfg.time_index, it);
    end

    if isempty(sigma_vec)
        sigma_vec = sigma(:);
        x_vec = x(:, 1);
        t_selected = t(cfg.time_index);
    end

    sigma_idx = find(sigma_vec == max(sigma_vec), 1, 'first');
    eta_phases(idx, :) = squeeze(eta(cfg.time_index, :, 1)); %#ok<SAGROW>
    phi_surface_phases(idx, :) = squeeze(phi(cfg.time_index, sigma_idx, :, 1)); %#ok<SAGROW>
    phit_corr_surface_phases(idx, :) = squeeze(phit_m(cfg.time_index, sigma_idx, :, 1)); %#ok<SAGROW>

    phit_uncorr = compute_uncorrected_phit_surface(phi(:, sigma_idx, :, 1), t);
    phit_uncorr_surface_phases(idx, :) = phit_uncorr(cfg.time_index, :); %#ok<SAGROW>
end

coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

eta_h = reconstruct_harmonics_1d_local(eta_phases, coef);
phi_h = reconstruct_harmonics_1d_local(phi_surface_phases, coef);
phit_corr_h = reconstruct_harmonics_1d_local(phit_corr_surface_phases, coef);
phit_uncorr_h = reconstruct_harmonics_1d_local(phit_uncorr_surface_phases, coef);

kp = 2 * pi / cfg.lambda;
if cfg.apply_x_filter
    eta_h = filter_harmonics_x_only_local(eta_h, x_vec, kp);
    phi_h = filter_harmonics_x_only_local(phi_h, x_vec, kp);
    phit_corr_h = filter_harmonics_x_only_local(phit_corr_h, x_vec, kp);
    phit_uncorr_h = filter_harmonics_x_only_local(phit_uncorr_h, x_vec, kp);
end

eta1 = eta_h(1, :).';
phi1 = phi_h(1, :).';
phit1_corr = phit_corr_h(1, :).';
phit1_uncorr = phit_uncorr_h(1, :).';
g_eta1 = -cfg.gravity * eta1;

stats = struct();
stats.corr_phitcorr_geta = corr_safe(phit1_corr, g_eta1);
stats.corr_phituncorr_geta = corr_safe(phit1_uncorr, g_eta1);
stats.scale_phitcorr_geta = least_squares_scale(eta1, phit1_corr);
stats.scale_phituncorr_geta = least_squares_scale(eta1, phit1_uncorr);

fprintf('\nPHIT DIAGNOSTICS AT SURFACE, t index = %d, t = %.6f s\n', cfg.time_index, t_selected);
fprintf('corr(phit_corr^(1), -g eta^(1))     = %.6f\n', stats.corr_phitcorr_geta);
fprintf('corr(phit_uncorr^(1), -g eta^(1))   = %.6f\n', stats.corr_phituncorr_geta);
fprintf('LS scale: phit_corr^(1) ~ c * eta^(1), c = %.6f\n', stats.scale_phitcorr_geta);
fprintf('LS scale: phit_uncorr^(1) ~ c * eta^(1), c = %.6f\n', stats.scale_phituncorr_geta);

x_plot = (x_vec - 0.5 * (x_vec(1) + x_vec(end))) / cfg.lambda;
fig = figure('Color', 'w', 'Position', [120 100 1450 920]);
tiledlayout(fig, 3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(x_plot, phit1_corr, 'k-', 'LineWidth', 1.8); hold on;
plot(x_plot, g_eta1, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8);
grid on; box on;
title('(a) Corrected \phi_t^{(1)} vs -g\eta^{(1)}');
legend({'Corrected \phi_t^{(1)}', '-g\eta^{(1)}'}, 'Location', 'best');
xlabel('x / \lambda');
ylabel('m^2/s^2');

nexttile;
plot(x_plot, phit1_uncorr, 'k-', 'LineWidth', 1.8); hold on;
plot(x_plot, g_eta1, '--', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8);
grid on; box on;
title('(b) Uncorrected D_t\phi^{(1)} vs -g\eta^{(1)}');
legend({'Uncorrected D_t\phi^{(1)}', '-g\eta^{(1)}'}, 'Location', 'best');
xlabel('x / \lambda');
ylabel('m^2/s^2');

nexttile;
plot(x_plot, phit1_corr, 'k-', 'LineWidth', 1.8); hold on;
plot(x_plot, phit1_uncorr, '--', 'Color', [0.12 0.39 0.71], 'LineWidth', 1.8);
plot(x_plot, g_eta1, ':', 'Color', [0.82 0.24 0.14], 'LineWidth', 1.8);
grid on; box on;
title('(c) Corrected vs uncorrected \phi_t^{(1)} and -g\eta^{(1)}');
legend({'Corrected \phi_t^{(1)}', 'Uncorrected D_t\phi^{(1)}', '-g\eta^{(1)}'}, 'Location', 'best');
xlabel('x / \lambda');
ylabel('m^2/s^2');

exportgraphics(fig, fullfile(pwd, 'processed_boundkinematics', 'OW3D_phit_diagnostic.png'), 'Resolution', 300);

function phit_uncorr = compute_uncorrected_phit_surface(phi_surface, t_vec)
    nt = size(phi_surface, 1);
    phi_surface = squeeze(phi_surface);
    nx = size(phi_surface, 2);
    dt = mean(diff(t_vec(:)));
    if ~isfinite(dt) || dt <= 0
        error('Invalid time vector for phit diagnostic.');
    end
    alpha = 2;
    r = 2 * alpha + 1;
    c = build_stencil_even_local(alpha, 1);
    dt_matrix = zeros(nt, nt);
    for j = alpha + 1:nt - alpha
        dt_matrix(j, j - alpha:j + alpha) = c(:, alpha + 1).';
    end
    for j = 1:alpha
        dt_matrix(j, :) = 0;
        dt_matrix(j, 1:r) = c(:, j).';
        dt_matrix(nt - j + 1, :) = 0;
        dt_matrix(nt - j + 1, nt - r + 1:nt) = c(:, r - j + 1).';
    end
    dt_matrix = dt_matrix / dt;
    phit_uncorr = zeros(nt, nx);
    for ip = 1:nx
        phi_col = phi_surface(:, ip);
        phit_uncorr(:, ip) = dt_matrix * phi_col;
    end
end

function h = reconstruct_harmonics_1d_local(phase_data, coef)
    n_phase = size(phase_data, 1);
    nx = size(phase_data, 2);
    h = zeros(4, nx);
    for ix = 1:nx
        y = phase_data(:, ix);
        analytic_part = imag(hilbert(y));
        state = [y(1); analytic_part(1); y(2); analytic_part(2); y(3); analytic_part(3); y(4); analytic_part(4)];
        h(:, ix) = coef * state;
    end
end

function filtered = filter_harmonics_x_only_local(harmonics, x_vec, kp)
    filtered = zeros(size(harmonics));
    for n = 1:size(harmonics, 1)
        filtered(n, :) = frequency_filtering_1d_local(harmonics(n, :), x_vec, kp, n);
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

function c = corr_safe(a, b)
    a = a(:);
    b = b(:);
    if all(abs(a) < eps) || all(abs(b) < eps)
        c = NaN;
        return;
    end
    c = corr(a, b, 'Rows', 'complete');
end

function s = least_squares_scale(a, b)
    a = a(:);
    b = b(:);
    denom = a' * a;
    if denom <= eps
        s = NaN;
    else
        s = (a' * b) / denom;
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
