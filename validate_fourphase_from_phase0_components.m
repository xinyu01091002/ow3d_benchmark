% validate_fourphase_from_phase0_components.m
% Validate a cheaper four-phase strategy without any spectral filtering.
%
% This version uses MF12 directly:
%   - order = 1  -> first-order field
%   - order = 2  -> first + second-order field
%   - order = 3  -> first + second + third-order field
%
% Then the perturbation-order increments are recovered by subtraction:
%   eta_1 = eta(order=1)
%   eta_2 = eta(order=2) - eta(order=1)
%   eta_3 = eta(order=3) - eta(order=2)
%
% The questions being tested are:
%   1) Can phase = 90/180/270 be reproduced from phase = 0 by shifting
%      eta_1, eta_2, eta_3 with 1*phi, 2*phi, 3*phi respectively?
%   2) If eta_2 is kept as a direct MF12 computation, can eta_3 alone be
%      reused cheaply from phase = 0 via a 3*phi phase shift?
%
% This version also diagnoses the second-order increment by splitting eta_2
% into a low-wavenumber subharmonic part and a 2kp superharmonic part.

clc;
clear;
close all;

CFG = struct();

% -------------------- Case setup --------------------
CFG.g = 9.81;
CFG.kp = 0.0279;
CFG.kd = 2;
CFG.Akp = 0.12;
CFG.Alpha = 8.0;
CFG.heading_deg = 0;
CFG.spread_deg = 25;
CFG.energy_keep_frac = 0.9999;
CFG.max_components = 1600;

% -------------------- Domain / sampling --------------------
CFG.Lx_lambda = 5;
CFG.Ly_lambda = 5;
CFG.Nx = 513;
CFG.Ny = 128;
CFG.t_eval_periods = 0.0;
CFG.phase_test_deg = [0, 90, 180, 270];
CFG.plot_compare_deg = [90, 180, 270];
CFG.second_subharmonic_cut = 0.75; % keep |k| <= 0.75 kp
CFG.second_super_sigma = 0.22;     % Gaussian width around 2 kp, in units of kp

setup_mf12_paths();

h = CFG.kd / CFG.kp;
lambda_p = 2 * pi / CFG.kp;
Lx = CFG.Lx_lambda * lambda_p;
Ly = CFG.Ly_lambda * lambda_p;
t_eval = CFG.t_eval_periods * 2 * pi / sqrt(CFG.g * CFG.kp * tanh(CFG.kp * h));

[kx, ky, amp_base] = build_directional_group_spectrum(CFG, Lx, Ly);

xf = 0.5 * Lx;
yf = 0.5 * Ly;
omega_lin = sqrt(CFG.g * hypot(kx, ky) .* tanh(h * hypot(kx, ky)));
phase_focus_0 = -(kx * xf + ky * yf);
a0 = amp_base .* cos(phase_focus_0);
b0 = amp_base .* sin(phase_focus_0);

[eta0, phi0] = reconstruct_by_order(CFG.g, h, a0, b0, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

results = struct([]);

for i = 1:numel(CFG.phase_test_deg)
    phi_deg = CFG.phase_test_deg(i);
    phi_rad = deg2rad(phi_deg);

    phase_direct = phase_focus_0 + phi_rad;
    a_dir = amp_base .* cos(phase_direct);
    b_dir = amp_base .* sin(phase_direct);

    [eta_direct, phi_direct] = reconstruct_by_order( ...
        CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

    eta_syn_pos = synthesize_from_phase0(eta0, Lx, Ly, phi_rad, +1);
    eta_syn_neg = synthesize_from_phase0(eta0, Lx, Ly, phi_rad, -1);
    phi_syn_pos = synthesize_from_phase0(phi0, Lx, Ly, phi_rad, +1);
    phi_syn_neg = synthesize_from_phase0(phi0, Lx, Ly, phi_rad, -1);

    results(i).phi_deg = phi_deg;
    results(i).eta_total_pos_err = max(abs(eta_direct.total(:) - eta_syn_pos.total(:)));
    results(i).eta_total_neg_err = max(abs(eta_direct.total(:) - eta_syn_neg.total(:)));
    results(i).eta_1_pos_err = max(abs(eta_direct.eta1(:) - eta_syn_pos.eta1(:)));
    results(i).eta_1_neg_err = max(abs(eta_direct.eta1(:) - eta_syn_neg.eta1(:)));
    results(i).eta_2_pos_err = max(abs(eta_direct.eta2(:) - eta_syn_pos.eta2(:)));
    results(i).eta_2_neg_err = max(abs(eta_direct.eta2(:) - eta_syn_neg.eta2(:)));
    results(i).eta_3_pos_err = max(abs(eta_direct.eta3(:) - eta_syn_pos.eta3(:)));
    results(i).eta_3_neg_err = max(abs(eta_direct.eta3(:) - eta_syn_neg.eta3(:)));
    results(i).phi_total_pos_err = max(abs(phi_direct.total(:) - phi_syn_pos.total(:)));
    results(i).phi_total_neg_err = max(abs(phi_direct.total(:) - phi_syn_neg.total(:)));
end

fprintf('\nOrder-increment validation summary\n');
fprintf('phi(deg) | total(+phi) | total(-phi) | eta1(+/-)     | eta2(+/-)     | eta3(+/-)\n');
fprintf('----------------------------------------------------------------------------------\n');
for i = 1:numel(results)
    fprintf('%7d | %11.3e | %11.3e | %6.2e / %6.2e | %6.2e / %6.2e | %6.2e / %6.2e\n', ...
        results(i).phi_deg, ...
        results(i).eta_total_pos_err, ...
        results(i).eta_total_neg_err, ...
        results(i).eta_1_pos_err, results(i).eta_1_neg_err, ...
        results(i).eta_2_pos_err, results(i).eta_2_neg_err, ...
        results(i).eta_3_pos_err, results(i).eta_3_neg_err);
end

fprintf('\nHybrid strategy diagnostic (direct order<=2, shifted order-3 only)\n');
fprintf('phi(deg) | hybrid total(+phi) | hybrid total(-phi) | eta3(+phi) | eta3(-phi)\n');
fprintf('-------------------------------------------------------------------------------\n');
for i = 1:numel(CFG.phase_test_deg)
    phi_deg = CFG.phase_test_deg(i);
    phi_rad = deg2rad(phi_deg);

    phase_direct = phase_focus_0 + phi_rad;
    a_dir = amp_base .* cos(phase_direct);
    b_dir = amp_base .* sin(phase_direct);
    [eta_direct, phi_direct] = reconstruct_by_order( ...
        CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

    eta3_pos = shift_field_in_fft(eta0.eta3, Lx, Ly, phi_rad, 3, +1);
    eta3_neg = shift_field_in_fft(eta0.eta3, Lx, Ly, phi_rad, 3, -1);
    phi3_pos = shift_field_in_fft(phi0.eta3, Lx, Ly, phi_rad, 3, +1);
    phi3_neg = shift_field_in_fft(phi0.eta3, Lx, Ly, phi_rad, 3, -1);

    eta_hybrid_pos = eta_direct.eta1 + eta_direct.eta2 + eta3_pos;
    eta_hybrid_neg = eta_direct.eta1 + eta_direct.eta2 + eta3_neg;
    phi_hybrid_pos = phi_direct.eta1 + phi_direct.eta2 + phi3_pos;
    phi_hybrid_neg = phi_direct.eta1 + phi_direct.eta2 + phi3_neg;

    fprintf('%7d | %18.3e | %18.3e | %10.3e | %10.3e\n', ...
        phi_deg, ...
        max(abs(eta_direct.total(:) - eta_hybrid_pos(:))), ...
        max(abs(eta_direct.total(:) - eta_hybrid_neg(:))), ...
        max(abs(eta_direct.eta3(:) - eta3_pos(:))), ...
        max(abs(eta_direct.eta3(:) - eta3_neg(:))));

    results(i).hybrid_eta_total_pos_err = max(abs(eta_direct.total(:) - eta_hybrid_pos(:)));
    results(i).hybrid_eta_total_neg_err = max(abs(eta_direct.total(:) - eta_hybrid_neg(:)));
    results(i).hybrid_phi_total_pos_err = max(abs(phi_direct.total(:) - phi_hybrid_pos(:)));
    results(i).hybrid_phi_total_neg_err = max(abs(phi_direct.total(:) - phi_hybrid_neg(:)));
end

fprintf('\nSecond-order split diagnostic\n');
fprintf('phi(deg) | eta2_sub(0-shift) | eta2_sup(+phi) | eta2_sup(-phi)\n');
fprintf('---------------------------------------------------------------\n');
for i = 1:numel(CFG.phase_test_deg)
    phi_deg = CFG.phase_test_deg(i);
    phi_rad = deg2rad(phi_deg);

    phase_direct = phase_focus_0 + phi_rad;
    a_dir = amp_base .* cos(phase_direct);
    b_dir = amp_base .* sin(phase_direct);
    eta_direct = reconstruct_by_order(CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

    eta2_direct_split = split_second_order_components(eta_direct.eta2, Lx, Ly, CFG.kp, CFG);
    eta2_ref_split = split_second_order_components(eta0.eta2, Lx, Ly, CFG.kp, CFG);
    eta2_sup_pos = shift_field_in_fft(eta2_ref_split.sup, Lx, Ly, phi_rad, 2, +1);
    eta2_sup_neg = shift_field_in_fft(eta2_ref_split.sup, Lx, Ly, phi_rad, 2, -1);

    fprintf('%7d | %16.3e | %14.3e | %14.3e\n', ...
        phi_deg, ...
        max(abs(eta2_direct_split.sub(:) - eta2_ref_split.sub(:))), ...
        max(abs(eta2_direct_split.sup(:) - eta2_sup_pos(:))), ...
        max(abs(eta2_direct_split.sup(:) - eta2_sup_neg(:))));
end

x = (0:CFG.Nx-1) * (Lx / CFG.Nx);
y = (0:CFG.Ny-1) * (Ly / CFG.Ny);
[~, iyc] = min(abs(y - 0.5 * Ly));
compare_phases = intersect(CFG.plot_compare_deg, CFG.phase_test_deg, 'stable');

if ~isempty(compare_phases)
    fig = figure('Color', 'w', 'Position', [120 120 1500 900]);
    tiledlayout(2, numel(compare_phases), 'Padding', 'compact', 'TileSpacing', 'compact');

    for j = 1:numel(compare_phases)
        phi_deg = compare_phases(j);
        phi_rad = deg2rad(phi_deg);

        phase_direct = phase_focus_0 + phi_rad;
        a_dir = amp_base .* cos(phase_direct);
        b_dir = amp_base .* sin(phase_direct);
        eta_direct = reconstruct_by_order(CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);
        eta_syn_pos = synthesize_from_phase0(eta0, Lx, Ly, phi_rad, +1);
        eta2_direct_split = split_second_order_components(eta_direct.eta2, Lx, Ly, CFG.kp, CFG);
        eta2_ref_split = split_second_order_components(eta0.eta2, Lx, Ly, CFG.kp, CFG);
        eta2_sup_pos = shift_field_in_fft(eta2_ref_split.sup, Lx, Ly, phi_rad, 2, +1);

        nexttile(j);
        plot(x, eta_direct.eta2(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct from phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta_syn_pos.eta2(iyc, :), 'r--', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('From phase=0 shifted to %d^\\circ', phi_deg));
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Second-order increment centerline, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta_2 (m)');

        nexttile(numel(compare_phases) + j);
        plot(x, eta_direct.eta3(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct from phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta_syn_pos.eta3(iyc, :), 'r--', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('From phase=0 shifted to %d^\\circ', phi_deg));
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Third-order increment centerline, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta_3 (m)');
    end
end

if ~isempty(compare_phases)
    fig_hybrid = figure('Color', 'w', 'Position', [160 160 1500 900]);
    tiledlayout(2, numel(compare_phases), 'Padding', 'compact', 'TileSpacing', 'compact');

    for j = 1:numel(compare_phases)
        phi_deg = compare_phases(j);
        phi_rad = deg2rad(phi_deg);

        phase_direct = phase_focus_0 + phi_rad;
        a_dir = amp_base .* cos(phase_direct);
        b_dir = amp_base .* sin(phase_direct);
        [eta_direct, ~] = reconstruct_by_order(CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

        eta3_pos = shift_field_in_fft(eta0.eta3, Lx, Ly, phi_rad, 3, +1);
        eta_hybrid_pos = eta_direct.eta1 + eta_direct.eta2 + eta3_pos;

        nexttile(j);
        plot(x, eta_direct.eta3(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct third-order, phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta3_pos(iyc, :), 'r--', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('From phase=0, third-order shifted to %d^\\circ', phi_deg));
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Third-order only reuse, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta_3 (m)');

        nexttile(numel(compare_phases) + j);
        plot(x, eta_direct.total(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct total, phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta_hybrid_pos(iyc, :), 'r--', 'LineWidth', 1.4, ...
            'DisplayName', 'Hybrid: direct order<=2 + shifted order-3');
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Hybrid total field, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta (m)');
    end
end

if ~isempty(compare_phases)
    fig_eta2 = figure('Color', 'w', 'Position', [140 140 1500 900]);
    tiledlayout(2, numel(compare_phases), 'Padding', 'compact', 'TileSpacing', 'compact');

    eta2_ref_split = split_second_order_components(eta0.eta2, Lx, Ly, CFG.kp, CFG);

    for j = 1:numel(compare_phases)
        phi_deg = compare_phases(j);
        phi_rad = deg2rad(phi_deg);

        phase_direct = phase_focus_0 + phi_rad;
        a_dir = amp_base .* cos(phase_direct);
        b_dir = amp_base .* sin(phase_direct);
        eta_direct = reconstruct_by_order(CFG.g, h, a_dir, b_dir, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);
        eta2_direct_split = split_second_order_components(eta_direct.eta2, Lx, Ly, CFG.kp, CFG);
        eta2_sup_pos = shift_field_in_fft(eta2_ref_split.sup, Lx, Ly, phi_rad, 2, +1);

        nexttile(j);
        plot(x, eta2_direct_split.sup(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct super, phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta2_sup_pos(iyc, :), 'r--', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('From phase=0 super shifted to %d^\\circ', phi_deg));
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Second-order superharmonic, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta_{2,sup} (m)');

        nexttile(numel(compare_phases) + j);
        plot(x, eta2_direct_split.sub(iyc, :), 'k-', 'LineWidth', 1.4, ...
            'DisplayName', sprintf('Direct sub, phase=%d^\\circ', phi_deg));
        hold on;
        plot(x, eta2_ref_split.sub(iyc, :), 'g--', 'LineWidth', 1.4, ...
            'DisplayName', 'Phase=0 sub reused');
        hold off;
        grid on;
        box on;
        legend('Location', 'best');
        title(sprintf('Second-order subharmonic, %d^\\circ', phi_deg));
        xlabel('x (m)');
        ylabel('\eta_{2,sub} (m)');
    end
end

% -------------------- Local helpers --------------------
function setup_mf12_paths()
    source_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
    if ~isfolder(source_dir)
        error('MF12 source directory not found: %s', source_dir);
    end
    addpath(source_dir);
end

function [kx, ky, amp] = build_directional_group_spectrum(CFG, Lx, Ly)
    dkx = 2 * pi / Lx;
    dky = 2 * pi / Ly;
    kx_idx_all = (-floor(CFG.Nx / 2)):(ceil(CFG.Nx / 2) - 1);
    ky_idx_all = (-floor(CFG.Ny / 2)):(ceil(CFG.Ny / 2) - 1);
    [KXI, KYI] = meshgrid(kx_idx_all, ky_idx_all);

    kx_all = KXI(:) * dkx;
    ky_all = KYI(:) * dky;
    kmag_all = hypot(kx_all, ky_all);
    theta_all = atan2(ky_all, kx_all);

    keep = (kx_all > 0) | (kx_all == 0 & ky_all > 0);
    kx_all = kx_all(keep);
    ky_all = ky_all(keep);
    kmag_all = kmag_all(keep);
    theta_all = theta_all(keep);

    kw_right = sqrt(CFG.kp^2 / (2 * log(10^CFG.Alpha)));
    kw_vec = kw_right * ones(size(kmag_all));
    kw_vec(kmag_all <= CFG.kp) = 0.004606;
    Sk = exp(-((kmag_all - CFG.kp).^2) ./ (2 * kw_vec.^2));
    D = gaussian_spreading(theta_all - deg2rad(CFG.heading_deg), CFG.spread_deg);
    W = Sk .* D;

    valid = W > 1e-10 * max(W);
    kx_all = kx_all(valid);
    ky_all = ky_all(valid);
    W = W(valid);

    [W_sorted, idx_sort] = sort(W, 'descend');
    cum_energy = cumsum(W_sorted);
    n_keep = find(cum_energy >= CFG.energy_keep_frac * cum_energy(end), 1, 'first');
    n_keep = min(n_keep, CFG.max_components);
    idx = idx_sort(1:n_keep);

    kx = kx_all(idx).';
    ky = ky_all(idx).';
    amp = W(idx).';
    amp = amp * ((CFG.Akp / CFG.kp) / max(sum(amp), eps));
end

function D = gaussian_spreading(theta, spread_angle_deg)
    theta_wrapped = angle(exp(1i * theta));
    sigma = deg2rad(spread_angle_deg);
    D = exp(-0.5 * (theta_wrapped / max(sigma, eps)).^2);
end

function [eta_parts, phi_parts] = reconstruct_by_order(g, h, a, b, kx, ky, Lx, Ly, Nx, Ny, t_eval)
    coeffs1 = mf12_spectral_coefficients(1, g, h, a, b, kx, ky, 0, 0);
    coeffs2 = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, 0, 0);
    coeffs3 = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, 0, 0);

    [eta_o1, phi_o1] = mf12_spectral_surface(coeffs1, Lx, Ly, Nx, Ny, t_eval);
    [eta_o2, phi_o2] = mf12_spectral_surface(coeffs2, Lx, Ly, Nx, Ny, t_eval);
    [eta_o3, phi_o3] = mf12_spectral_surface(coeffs3, Lx, Ly, Nx, Ny, t_eval);

    eta_parts = struct();
    eta_parts.eta1 = eta_o1;
    eta_parts.eta2 = eta_o2 - eta_o1;
    eta_parts.eta3 = eta_o3 - eta_o2;
    eta_parts.total = eta_parts.eta1 + eta_parts.eta2 + eta_parts.eta3;

    phi_parts = struct();
    phi_parts.eta1 = phi_o1;
    phi_parts.eta2 = phi_o2 - phi_o1;
    phi_parts.eta3 = phi_o3 - phi_o2;
    phi_parts.total = phi_parts.eta1 + phi_parts.eta2 + phi_parts.eta3;
end

function parts_shifted = synthesize_from_phase0(parts0, Lx, Ly, phi_rad, sign_flag)
    parts_shifted = struct();
    parts_shifted.eta1 = shift_field_in_fft(parts0.eta1, Lx, Ly, phi_rad, 1, sign_flag);
    parts_shifted.eta2 = shift_field_in_fft(parts0.eta2, Lx, Ly, phi_rad, 2, sign_flag);
    parts_shifted.eta3 = shift_field_in_fft(parts0.eta3, Lx, Ly, phi_rad, 3, sign_flag);
    parts_shifted.total = parts_shifted.eta1 + parts_shifted.eta2 + parts_shifted.eta3;
end

function eta2_split = split_second_order_components(eta2, Lx, Ly, kp, CFG)
    eta2_split = struct();
    eta2_split.sub = lowpass_component(eta2, Lx, Ly, CFG.second_subharmonic_cut * kp);
    eta2_split.sup = radial_band_component(eta2, Lx, Ly, 2.0 * kp, CFG.second_super_sigma * kp);
    eta2_split.remainder = eta2 - eta2_split.sub - eta2_split.sup;
end

function comp = lowpass_component(field_in, Lx, Ly, k_cut)
    comp = masked_fft_component(field_in, Lx, Ly, @(K) double(K <= k_cut));
end

function comp = radial_band_component(field_in, Lx, Ly, k_target, sigma)
    comp = masked_fft_component(field_in, Lx, Ly, ...
        @(K) exp(-((K - k_target).^2) / (2 * sigma^2)));
end

function comp = masked_fft_component(field_in, Lx, Ly, mask_fn)
    [Ny, Nx] = size(field_in);
    kx_idx = (-floor(Nx / 2)):(ceil(Nx / 2) - 1);
    ky_idx = (-floor(Ny / 2)):(ceil(Ny / 2) - 1);
    dkx = 2 * pi / Lx;
    dky = 2 * pi / Ly;
    [KXI, KYI] = meshgrid(kx_idx, ky_idx);
    K = hypot(KXI * dkx, KYI * dky);

    F = fftshift(fft2(field_in));
    mask = mask_fn(K);
    comp = ifft2(ifftshift(F .* mask));
    if isreal(field_in)
        comp = real(comp);
    end
end

function field_shifted = shift_field_in_fft(field0, Lx, Ly, phase_shift, harmonic_order, sign_flag)
    [Ny, Nx] = size(field0);
    kx_idx = (-floor(Nx / 2)):(ceil(Nx / 2) - 1);
    ky_idx = (-floor(Ny / 2)):(ceil(Ny / 2) - 1);
    dkx = 2 * pi / Lx;
    dky = 2 * pi / Ly;
    [KXI, KYI] = meshgrid(kx_idx, ky_idx);
    kx = KXI * dkx;
    ky = KYI * dky;

    pos_mask = (kx > 0) | (kx == 0 & ky > 0);
    neg_mask = (kx < 0) | (kx == 0 & ky < 0);

    F = fftshift(fft2(field0));
    F(pos_mask) = F(pos_mask) .* exp(1i * sign_flag * harmonic_order * phase_shift);
    F(neg_mask) = F(neg_mask) .* exp(-1i * sign_flag * harmonic_order * phase_shift);

    field_shifted = ifft2(ifftshift(F));
    if isreal(field0)
        field_shifted = real(field_shifted);
    end
end
