%% Diagnose MF12 first-order error against OW3D first harmonic at z = 0

clc; clear;

current_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(current_dir);
addpath(fullfile(project_dir, 'irregularWavesMF12', 'Source'));

raw_path = fullfile(project_dir, 'processed_boundkinematics', 'OW3D_boundkinematics_raw_tidx_0066.mat');
S = load(raw_path);

phase_names = sort(fieldnames(S.raw_snapshot));
x_vec = S.raw_meta.x(:);
sigma_vec = S.raw_meta.sigma(:);
[~, sigma_idx] = max(sigma_vec);
h = mean(S.raw_snapshot.(phase_names{1}).h(:), 'omitnan');
g = 9.81;

eta_phases = zeros(numel(phase_names), numel(x_vec));
u_phases = zeros(numel(phase_names), numel(x_vec));
w_phases = zeros(numel(phase_names), numel(x_vec));

for idx = 1:numel(phase_names)
    snapshot = S.raw_snapshot.(phase_names{idx});
    eta_phases(idx, :) = squeeze(snapshot.eta(1, :));
    u_phases(idx, :) = squeeze(snapshot.u(sigma_idx, :));
    w_phases(idx, :) = squeeze(snapshot.w(sigma_idx, :));
end

coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

eta_h = reconstruct_harmonics_1d_local(eta_phases, coef);
u_h = reconstruct_harmonics_1d_local(u_phases, coef);
w_h = reconstruct_harmonics_1d_local(w_phases, coef);

eta11 = eta_h(1, :).';
[mf12_eta1, mf12_u1, mf12_w1, mf12_meta] = compute_mf12_first_order_surface(eta11, x_vec, h, g);

fprintf('MF12 first-order reconstruction from eta11 snapshot\n');
fprintf('Retained components = %d\n', mf12_meta.num_components);
fprintf('eta1 corr = %.6f, rmse = %.6e, peak ratio = %.6f\n', metric_corr(eta_h(1, :).', mf12_eta1), ...
    metric_rmse(eta_h(1, :).', mf12_eta1), metric_peak_ratio(eta_h(1, :).', mf12_eta1));
fprintf('u1 corr   = %.6f, rmse = %.6e, peak ratio = %.6f\n', metric_corr(u_h(1, :).', mf12_u1), ...
    metric_rmse(u_h(1, :).', mf12_u1), metric_peak_ratio(u_h(1, :).', mf12_u1));
fprintf('w1 corr   = %.6f, rmse = %.6e, peak ratio = %.6f\n', metric_corr(w_h(1, :).', mf12_w1), ...
    metric_rmse(w_h(1, :).', mf12_w1), metric_peak_ratio(w_h(1, :).', mf12_w1));

function h = reconstruct_harmonics_1d_local(phase_data, coef)
    analytic_part = hilbert(phase_data.').';
    all_fields = zeros(8, size(phase_data, 2));
    all_fields(1:2:end, :) = real(phase_data);
    all_fields(2:2:end, :) = -imag(analytic_part);
    h = zeros(4, size(phase_data, 2));
    for n = 1:4
        h(n, :) = coef(n, :) * all_fields;
    end
end

function [eta1, u1, w1, meta] = compute_mf12_first_order_surface(eta11, x_vec, h, g)
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
    kx = kx_grid(positive_idx).';
    ky = zeros(size(kx));
    a = 2 * real(fft_eta(positive_idx)).';
    b = 2 * imag(fft_eta(positive_idx)).';

    coeffs = mf12_direct_coefficients(1, g, h, a, b, kx, ky, 0, 0, 0);
    [u, ~, w, ~, phi] = kinematicsMF12(1, coeffs, x_vec.', 0, 0, 0); %#ok<ASGLU>

    theta = coeffs.omega(:) .* 0 - coeffs.kx(:) .* x_vec.' - coeffs.ky(:) .* 0;
    eta_matrix = coeffs.a(:) .* cos(theta) + coeffs.b(:) .* sin(theta);
    eta1 = sum(eta_matrix, 1).';

    u1 = u(:);
    w1 = w(:);
    meta = struct('num_components', numel(kx), 'coeffs', coeffs, 'phi', phi);
end

function out = metric_corr(a, b)
    a = a(:);
    b = b(:);
    cc = corrcoef(a, b);
    out = cc(1, 2);
end

function out = metric_rmse(a, b)
    a = a(:);
    b = b(:);
    out = sqrt(mean((a - b).^2));
end

function out = metric_peak_ratio(ref, candidate)
    ref = ref(:);
    candidate = candidate(:);
    out = max(abs(candidate)) / max(abs(ref));
end
