%% Diagnose MF12-native second subharmonic metrics at z = 0

clc; clear;

current_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(current_dir);
addpath(fullfile(project_dir, 'irregularWavesMF12', 'Source'));

S = load(fullfile(project_dir, 'processed_boundkinematics', 'OW3D_boundkinematics_raw_tidx_0066.mat'));
phase_names = sort(fieldnames(S.raw_snapshot));
x = S.raw_meta.x(:);
sigma_vec = S.raw_meta.sigma(:);
[~, sigma_idx] = max(sigma_vec);

eta_phases = zeros(numel(phase_names), numel(x));
u_phases = zeros(numel(phase_names), numel(x));
w_phases = zeros(numel(phase_names), numel(x));

for idx = 1:numel(phase_names)
    snapshot = S.raw_snapshot.(phase_names{idx});
    eta_phases(idx, :) = squeeze(snapshot.eta(1, :));
    u_phases(idx, :) = squeeze(snapshot.u(sigma_idx, :));
    w_phases(idx, :) = squeeze(snapshot.w(sigma_idx, :));
end

coef = [ ...
    0.25  0    -0.25  0     0    -0.25  0     0.25; ...
    0.25 -0.25  0.25 -0.25  0     0      0     0; ...
    0.25  0    -0.25  0     0     0.25  0    -0.25; ...
    0.25  0.25  0.25  0.25  0     0      0     0];
eta_h = reconstruct_harmonics_1d_local(eta_phases, coef);

kp = 2 * pi / 225;
h = mean(S.raw_snapshot.(phase_names{1}).h(:), 'omitnan');
g = 9.81;

ow3d_u20 = lowpass_wavenumber_component_local(mean(u_phases, 1).', x, 1.0 * kp, 0.35 * kp);
ow3d_w20 = lowpass_wavenumber_component_local(mean(w_phases, 1).', x, 1.0 * kp, 0.35 * kp);
eta_linear = eta_h(1, :).';
mf12 = compute_mf12_second_order_filtered_surface(eta_linear, x, h, g, 1.0 * kp, 0.35 * kp, 0.99);
fprintf('Retained %d linear components (energy keep = %.2f)\n', numel(mf12.linear_indices), mf12.linear_energy_keep);

fprintf('MF12 max eta20 = %.6e\n', max(abs(mf12.eta20)));
fprintf('MF12 max phi20 = %.6e\n', max(abs(mf12.phi20)));
fprintf('MF12 max u20   = %.6e\n', max(abs(mf12.u20)));
fprintf('MF12 max w20   = %.6e\n', max(abs(mf12.w20)));

report_metrics('u20', ow3d_u20, mf12.u20);
report_metrics('w20', ow3d_w20, mf12.w20);

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

    theta2 = 2 * coeffs2.omega(:) .* 0 - coeffs2.kx_2(:) .* x_vec.' - coeffs2.ky_2(:) .* 0;
    eta2_self = sum(coeffs2.G_2(:) .* (coeffs2.A_2(:) .* cos(theta2) + coeffs2.B_2(:) .* sin(theta2)), 1).';

    eta2_pair = zeros(nx, 1);
    cnm = 0;
    for n = 1:coeffs2.N
        for m = (n + 1):coeffs2.N
            for pm = [1 -1] %#ok<NASGU>
                cnm = cnm + 1;
                theta_npm = coeffs2.omega_npm(cnm) .* 0 - coeffs2.kx_npm(cnm) .* x_vec - coeffs2.ky_npm(cnm) .* 0;
                eta2_pair = eta2_pair + coeffs2.G_npm(cnm) .* ...
                    (coeffs2.A_npm(cnm) .* cos(theta_npm) + coeffs2.B_npm(cnm) .* sin(theta_npm));
            end
        end
    end

    eta20 = eta2_self + eta2_pair;
    phi20 = phi2(:) - phi1(:);
    u20 = u2(:) - u1(:);
    w20 = w2(:) - w1(:);

    eta20 = lowpass_wavenumber_component_local(eta20, x_vec, k_cutoff, transition);
    phi20 = lowpass_wavenumber_component_local(phi20, x_vec, k_cutoff, transition);
    u20 = lowpass_wavenumber_component_local(u20, x_vec, k_cutoff, transition);
    w20 = lowpass_wavenumber_component_local(w20, x_vec, k_cutoff, transition);

    out = struct('eta20', eta20(:), 'phi20', phi20(:), 'u20', u20(:), 'w20', w20(:), ...
        'linear_indices', positive_idx(:), 'linear_energy_keep', energy_keep);
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

function report_metrics(name, reference, candidate)
    reference = reference(:);
    candidate = candidate(:);
    cc = corrcoef(reference, candidate);
    fprintf('%s corr       = %.6f\n', name, cc(1, 2));
    fprintf('%s rmse       = %.6e\n', name, sqrt(mean((reference - candidate).^2)));
    fprintf('%s peak ratio = %.6f\n', name, max(abs(candidate)) / max(abs(reference)));
end

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
