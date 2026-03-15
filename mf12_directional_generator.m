% mf12_directional_generator.m
% Generate third-order directional wave-group initial conditions for OceanWave3D
% using the MF12 spectral reconstruction workflow.

clc;
clear;
close all;

CFG = struct();

% -------------------- Sea state --------------------
CFG.g = 9.81;
CFG.kp = 0.00279;
CFG.Akp_list = [0.12];
CFG.Alpha_list = [8];
CFG.kd_list = [1];
CFG.phases_deg = 0:90:270;
CFG.single_case_only = true;

% -------------------- Directional wave-group definition --------------------
CFG.heading_deg = 0;
CFG.spread_deg = 15;
CFG.energy_keep_frac = 0.999;
CFG.max_components = 800;

% -------------------- Domain / timing --------------------
CFG.Tp = 12;
CFG.t_init_periods_list = [-40, -30, -20]; % Initial-condition times relative to focus in units of Tp.
CFG.t_end_periods = 5; % Target end time relative to focus in units of Tp.
CFG.t_focus = 0.0;
CFG.steps_per_period = 30;
CFG.Lx_lambda = 40; % Physical domain size in x: Lx = Lx_lambda * lambda_p. Increase this if you want a physically larger x-domain.
CFG.Ly_lambda = 15; % Physical domain size in y: Ly = Ly_lambda * lambda_p. Increase this if you want a physically larger y-domain.
CFG.Nx = 1025; % Spatial resolution in x only. Increasing Nx decreases dx = Lx/Nx, but does not change the physical domain size.
CFG.Ny = 257; % Spatial resolution in y only. Increasing Ny decreases dy = Ly/Ny, but does not change the physical domain size.
CFG.focus_edge_padding_fraction = 0.05; % Reserve this fraction of Lx at each x-boundary when auto-placing the focus point.

% -------------------- Output --------------------
CFG.output_dir = fullfile('directional initial condition', 'error_wave_separation');
CFG.store_surface_stride = 4; % Section 8 surface output only. 1 saves every time step, 2 every second step, etc. Negative values request ascii output.
CFG.surface_format = 1; % Keep the existing OW3D surface-output format used in this project.
CFG.reuse_third_order_from_phase0 = true; % Hybrid workflow: compute order<=2 directly for each phase, and reuse the phase=0 third-order increment via a 3*phi shift.
CFG.batch_purpose = 'Generate directional OW3D cases to measure how quickly the error wave separates from the main wave group.';
CFG.batch_notes = [ ...
    "Each case starts from an early MF12 initialization time and is marched to a common end time of +5Tp relative to focus."; ...
    "The t_init sweep is intended to show when the error wave becomes clearly separated from the dominant wave group."; ...
    "Each subdirectory is a runnable OW3D case and is named so different start times remain easy to compare during postprocessing."];

setup_mf12_paths();

lambda_p = 2 * pi / CFG.kp;
Lx = CFG.Lx_lambda * lambda_p;
Ly = CFG.Ly_lambda * lambda_p;
dx = Lx / CFG.Nx;
dy = Ly / CFG.Ny;
fprintf('MF12 directional generator\n');
fprintf('Domain: Lx=%.3f m, Ly=%.3f m, Nx=%d, Ny=%d\n', Lx, Ly, CFG.Nx, CFG.Ny);
fprintf('Grid: dx=%.3f m, dy=%.3f m\n', dx, dy);
fprintf('Initial condition times: [%s] Tp relative to focus\n', strjoin(string(CFG.t_init_periods_list), ', '));
fprintf('Target end time: %.2f Tp relative to focus\n', CFG.t_end_periods);

batch_root = fullfile(pwd, CFG.output_dir);
ensure_dir(batch_root);
write_batch_readme(fullfile(batch_root, 'readme'), CFG, Lx, Ly, dx, dy);

kd_list = CFG.kd_list;
Akp_list = CFG.Akp_list;
Alpha_list = CFG.Alpha_list;
phase_list = CFG.phases_deg;
t_init_list = CFG.t_init_periods_list;

if CFG.single_case_only
    kd_list = kd_list(1);
    Akp_list = Akp_list(1);
    Alpha_list = Alpha_list(1);
    phase_list = phase_list(1);
end

for kd = kd_list
    h = kd / CFG.kp;
    Tp_kd = CFG.Tp;
    dt_kd = Tp_kd / CFG.steps_per_period;

    for Akp = Akp_list
        for Alpha = Alpha_list
            [kx, ky, amp_base, meta_base] = build_directional_group_spectrum(CFG, Lx, Ly, h, Akp, Alpha);

            fprintf('kd=%.2f, Akp=%.3f, Alpha=%.1f: retained %d components\n', ...
                kd, Akp, Alpha, numel(kx));

            omega_lin = sqrt(CFG.g * hypot(kx, ky) .* tanh(h * hypot(kx, ky)));
            for t_init_periods = t_init_list
                duration_periods = CFG.t_end_periods - t_init_periods;
                if duration_periods <= 0
                    error('t_end_periods must be larger than t_init_periods. Got %.2f and %.2f.', ...
                        CFG.t_end_periods, t_init_periods);
                end

                t_eval = t_init_periods * Tp_kd;
                N_steps_kd = round(duration_periods * CFG.steps_per_period);
                [xf, yf, focus_x_fraction] = resolve_focus_point(CFG, Lx, Ly, t_init_periods);
                phase_focus_0 = -(kx * xf + ky * yf) + omega_lin * CFG.t_focus;
                meta = meta_base;
                meta.t_init_periods = t_init_periods;
                meta.t_end_periods = CFG.t_end_periods;
                meta.duration_periods = duration_periods;
                meta.focus_x = xf;
                meta.focus_y = yf;
                meta.focus_x_fraction = focus_x_fraction;

                phase0_cache = struct();
                if CFG.reuse_third_order_from_phase0
                    a0 = amp_base .* cos(phase_focus_0);
                    b0 = amp_base .* sin(phase_focus_0);
                    [eta0_parts, phi0_parts] = reconstruct_order_parts( ...
                        CFG.g, h, a0, b0, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);
                    phase0_cache.eta3 = eta0_parts.order3_increment;
                    phase0_cache.phi3 = phi0_parts.order3_increment;
                end

                for phi_shift_deg = phase_list
                    phase_focus = phase_focus_0 + deg2rad(phi_shift_deg);
                    a = amp_base .* cos(phase_focus);
                    b = amp_base .* sin(phase_focus);

                    [eta_parts, phi_parts] = reconstruct_order_parts( ...
                        CFG.g, h, a, b, kx, ky, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

                    eta_lin = eta_parts.order1_total;
                    phi_lin = phi_parts.order1_total;
                    eta_mf12 = eta_parts.order2_total + eta_parts.order3_increment;
                    phi_mf12 = phi_parts.order2_total + phi_parts.order3_increment;

                    if CFG.reuse_third_order_from_phase0 && phi_shift_deg ~= 0
                        phi_rad = deg2rad(phi_shift_deg);
                        eta_mf12 = eta_parts.order2_total + ...
                            shift_field_in_fft(phase0_cache.eta3, Lx, Ly, phi_rad, 3, +1);
                        phi_mf12 = phi_parts.order2_total + ...
                            shift_field_in_fft(phase0_cache.phi3, Lx, Ly, phi_rad, 3, +1);
                    end

                    % mf12_spectral_surface returns arrays of size [Ny x Nx].
                    % Transpose to match the historical OW3D export convention used in this repo: [Nx x Ny].
                    eta_out = eta_mf12.';
                    phi_out = phi_mf12.';

                    write_path = fullfile(pwd, CFG.output_dir, ...
                        sprintf('T_init%d_Tend%d_Tp_kd%.1f_spread_%d_heading_%d_Akp_%03d_alpha_%.1f_phi_%d', ...
                        round(t_init_periods), round(CFG.t_end_periods), kd, round(CFG.spread_deg), round(CFG.heading_deg), ...
                        round(Akp * 100), Alpha, phi_shift_deg));

                    ensure_dir(write_path);

                    write_ow3d_init(fullfile(write_path, 'OceanWave3D.init'), ...
                        eta_out, phi_out, dx, dy, Lx, Ly, dt_kd, Akp, phi_shift_deg);

                    write_ow3d_inp(fullfile(write_path, 'OceanWave3D.inp'), ...
                        CFG, Lx, Ly, h, CFG.Nx, CFG.Ny, N_steps_kd, dt_kd);

                    write_readme(fullfile(write_path, 'OW_readme.txt'), ...
                        CFG, meta, Akp, Alpha, kd, h, Tp_kd, dx, dy, t_eval, phi_shift_deg);

                    save_visualizations(write_path, eta_lin, eta_mf12, phi_lin, phi_mf12, Lx, Ly);

                    fprintf('  wrote: %s\n', write_path);
                end
            end
        end
    end
end

function setup_mf12_paths()
    source_dir = fullfile(fileparts(mfilename('fullpath')), 'irregularWavesMF12', 'Source');
    if ~isfolder(source_dir)
        error('MF12 source directory not found: %s', source_dir);
    end

    addpath(source_dir);
end

function [xf, yf, focus_x_fraction] = resolve_focus_point(CFG, Lx, Ly, t_init_periods)
    duration_periods = CFG.t_end_periods - t_init_periods;
    focus_x_fraction = -t_init_periods / duration_periods;

    pad = CFG.focus_edge_padding_fraction;
    focus_x_fraction = pad + (1 - 2 * pad) * focus_x_fraction;
    focus_x_fraction = min(max(focus_x_fraction, pad), 1 - pad);

    xf = focus_x_fraction * Lx;
    yf = 0.5 * Ly;
end

function [kx, ky, amp, meta] = build_directional_group_spectrum(CFG, Lx, Ly, h, Akp, Alpha)
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

    kw_left = 0.004606;
    kw_right = sqrt(CFG.kp^2 / (2 * log(10^Alpha)));
    kw_vec = kw_right * ones(size(kmag_all));
    kw_vec(kmag_all <= CFG.kp) = kw_left;
    Sk = exp(-((kmag_all - CFG.kp).^2) ./ (2 * kw_vec.^2));

    D = gaussian_spreading(theta_all - deg2rad(CFG.heading_deg), CFG.spread_deg);
    W = Sk .* D;

    valid = W > 1e-10 * max(W);
    kx_all = kx_all(valid);
    ky_all = ky_all(valid);
    kmag_all = kmag_all(valid);
    W = W(valid);

    [W_sorted, idx_sort] = sort(W, 'descend');
    cum_energy = cumsum(W_sorted);
    N_energy = find(cum_energy >= CFG.energy_keep_frac * cum_energy(end), 1, 'first');
    N_keep = min(N_energy, CFG.max_components);
    idx = idx_sort(1:N_keep);

    kx = kx_all(idx).';
    ky = ky_all(idx).';
    kmag = kmag_all(idx);
    amp = W(idx).';
    amp = amp * ((Akp / CFG.kp) / max(sum(amp), eps));

    meta = struct();
    meta.n_components = numel(kx);
    meta.energy_keep_frac = CFG.energy_keep_frac;
    meta.kmin = min(kmag);
    meta.kmax = max(kmag);
    meta.heading_deg = CFG.heading_deg;
    meta.spread_deg = CFG.spread_deg;
    meta.depth = h;
end

function D = gaussian_spreading(theta, spread_angle_deg)
    theta_wrapped = angle(exp(1i * theta));
    sigma = deg2rad(spread_angle_deg);
    D = exp(-0.5 * (theta_wrapped / max(sigma, eps)).^2);
end

function ensure_dir(path_str)
    if ~exist(path_str, 'dir')
        mkdir(path_str);
    end
end

function write_ow3d_init(file_name, eta, phi, dx, dy, Lx, Ly, dt, Akp, phi_shift_deg)
    nx = size(eta, 1);
    ny = size(eta, 2);

    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, ' H=%f nx=%d ny=%d dx=%f dy=%f akp=%f shift=%f', ...
        max(eta(:)), nx, ny, dx, dy, Akp, phi_shift_deg);
    fprintf(f, '\n%12e %12e %d %d %12e', Lx, Ly, nx, ny, dt);

    for ry = 1:ny
        for rx = 1:nx
            fprintf(f, '\n%12e %12e', eta(rx, ry), phi(rx, ry));
        end
    end
end

function write_ow3d_inp(file_name, CFG, Lx, Ly, h, nx, ny, n_steps, dt)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'Data for MF12 directional wave-group initialization %s <-\n', datestr(now, 0));
    fprintf(f, '-1 2 <-\n');
    fprintf(f, '%d %d %d %d %d %d 0 0 1 1 1 1 <-\n', round(Lx), round(Ly), round(h), nx, ny, 17);
    fprintf(f, '4 4 4 1 1 1 <-\n');
    fprintf(f, '%d %f 1 0. 1 <-\n', n_steps + 1, dt);
    fprintf(f, '9.81 <-\n');
    fprintf(f, '1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n');
    fprintf(f, '0.05 1.00 1.84 2 0 0 1 6 32 <-\n');
    fprintf(f, '%d %d <-\n', CFG.store_surface_stride, CFG.surface_format);
    fprintf(f, '1 0 <-\n');
    fprintf(f, '0 6 10 0.08 0.08 0.4 <-\n');
    fprintf(f, '0 8. 3 X 0.0 <-\n');
    fprintf(f, '0 0 <-\n');
    fprintf(f, '0 2.0 2 0 0 1 0 <-\n');
    fprintf(f, '0 <-\n');
    fprintf(f, '33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <-\n');
end

function write_readme(file_name, CFG, meta, Akp, Alpha, kd, h, Tp, dx, dy, t_eval, phi_shift_deg)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'MF12 third-order directional wave-group initialization\n');
    fprintf(f, 'Generated: %s\n', datestr(now, 0));
    fprintf(f, 'Akp=%.4f, Alpha=%.1f, kp=%.5f, kd=%.2f, h=%.4f\n', Akp, Alpha, CFG.kp, kd, h);
    fprintf(f, 'Tp=%.6f s, t_init=%.6f s, phase shift=%d deg\n', Tp, t_eval, phi_shift_deg);
    fprintf(f, 'Heading=%g deg, spreading sigma=%g deg\n', meta.heading_deg, meta.spread_deg);
    fprintf(f, 'Components=%d, energy keep frac=%.5f\n', meta.n_components, meta.energy_keep_frac);
    fprintf(f, 'k-range=[%.6f, %.6f] 1/m\n', meta.kmin, meta.kmax);
    fprintf(f, 'Focus point: x=%.6f m (%.4f Lx), y=%.6f m\n', meta.focus_x, meta.focus_x_fraction, meta.focus_y);
    fprintf(f, 'Grid: Nx=%d, Ny=%d, dx=%.6f m, dy=%.6f m\n', CFG.Nx, CFG.Ny, dx, dy);
    if CFG.reuse_third_order_from_phase0
        fprintf(f, 'Model: MF12 hybrid workflow (direct order<=2, phase=0 reused for order-3 increment)\n');
    else
        fprintf(f, 'Model: MF12 spectral coefficients/order-3, single directional wave group\n');
    end
    fprintf(f, 'Initial condition time relative to focus: %.2f Tp\n', meta.t_init_periods);
    fprintf(f, 'Target end time relative to focus: %.2f Tp\n', meta.t_end_periods);
    fprintf(f, 'Total OW3D duration after initialization: %.2f Tp\n', meta.duration_periods);
    fprintf(f, 'Surface output stride: every %d time step(s)\n', abs(CFG.store_surface_stride));
    fprintf(f, 'Kinematic output: disabled\n');
end

function write_batch_readme(file_name, CFG, Lx, Ly, dx, dy)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'Directional OW3D initial-condition batch\n');
    fprintf(f, 'Updated: %s\n\n', datestr(now, 0));

    fprintf(f, 'Purpose\n');
    fprintf(f, '%s\n\n', CFG.batch_purpose);

    fprintf(f, 'What this batch is trying to verify\n');
    for i = 1:numel(CFG.batch_notes)
        fprintf(f, '- %s\n', CFG.batch_notes(i));
    end
    fprintf(f, '\n');

    fprintf(f, 'Current workflow\n');
    if CFG.reuse_third_order_from_phase0
        fprintf(f, '- Direct MF12 reconstruction up to order 2 for each requested phase.\n');
        fprintf(f, '- Reuse the phase=0 third-order increment by applying a 3phi shift for the non-zero phase cases.\n');
    else
        fprintf(f, '- Direct MF12 order-3 reconstruction for every requested phase.\n');
    end
    fprintf(f, '- Surface-only OW3D output; no kinematic export.\n');
    fprintf(f, '- One subdirectory per generated case.\n\n');

    fprintf(f, 'Current settings snapshot\n');
    fprintf(f, '- output directory = %s\n', CFG.output_dir);
    fprintf(f, '- phases = [%s] deg\n', strjoin(string(CFG.phases_deg), ', '));
    fprintf(f, '- heading = %.1f deg\n', CFG.heading_deg);
    fprintf(f, '- spread = %.1f deg\n', CFG.spread_deg);
    fprintf(f, '- t_init sweep = [%s] Tp\n', strjoin(string(CFG.t_init_periods_list), ', '));
    fprintf(f, '- target end time = %.2f Tp\n', CFG.t_end_periods);
    fprintf(f, '- focus x-padding fraction = %.3f per boundary\n', CFG.focus_edge_padding_fraction);
    fprintf(f, '- domain = [%.3f, %.3f] m\n', Lx, Ly);
    fprintf(f, '- grid = [%d, %d], dx = %.3f m, dy = %.3f m\n', CFG.Nx, CFG.Ny, dx, dy);
    fprintf(f, '- surface stride = %d\n', CFG.store_surface_stride);
end

function save_visualizations(write_path, eta_lin, eta_total, phi_lin, phi_total, Lx, Ly)
    x = linspace(0, Lx, size(eta_lin, 2));
    y = linspace(0, Ly, size(eta_lin, 1));
    [~, iyc] = min(abs(y - 0.5 * Ly));

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1400 900]);
    tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    imagesc(x, y, eta_lin);
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('x (m)');
    ylabel('y (m)');
    title('Linear wave group: \eta');
    colorbar;

    nexttile;
    imagesc(x, y, eta_total);
    axis image;
    set(gca, 'YDir', 'normal');
    xlabel('x (m)');
    ylabel('y (m)');
    title('Third-order MF12: \eta');
    colorbar;

    nexttile;
    plot(x, eta_lin(iyc, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x, eta_total(iyc, :), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('x (m)');
    ylabel('\eta (m)');
    title('Centerline elevation');
    legend({'Linear', 'MF12 third-order'}, 'Location', 'best');

    nexttile;
    plot(x, phi_lin(iyc, :), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(x, phi_total(iyc, :), 'r--', 'LineWidth', 1.5);
    grid on;
    xlabel('x (m)');
    ylabel('\phi_s (m^2/s)');
    title('Centerline free-surface potential');
    legend({'Linear', 'MF12 third-order'}, 'Location', 'best');

    exportgraphics(fig, fullfile(write_path, 'mf12_directional_overview.png'), 'Resolution', 200);
    close(fig);

    fig_lin = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);
    surf(x, y, eta_lin, 'EdgeColor', 'none');
    view(35, 28);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('\eta (m)');
    title('Linear directional wave group');
    colorbar;
    camlight headlight;
    lighting gouraud;
    exportgraphics(fig_lin, fullfile(write_path, 'linear_wave_group_surface.png'), 'Resolution', 220);
    close(fig_lin);
end

function [eta_parts, phi_parts] = reconstruct_order_parts(g, h, a, b, kx, ky, Lx, Ly, Nx, Ny, t_eval)
    coeffs1 = mf12_spectral_coefficients(1, g, h, a, b, kx, ky, 0, 0);
    coeffs2 = mf12_spectral_coefficients(2, g, h, a, b, kx, ky, 0, 0);
    coeffs3 = mf12_spectral_coefficients(3, g, h, a, b, kx, ky, 0, 0);

    [eta1, phi1] = mf12_spectral_surface(coeffs1, Lx, Ly, Nx, Ny, t_eval);
    [eta2, phi2] = mf12_spectral_surface(coeffs2, Lx, Ly, Nx, Ny, t_eval);
    [eta3, phi3] = mf12_spectral_surface(coeffs3, Lx, Ly, Nx, Ny, t_eval);

    eta_parts = struct();
    eta_parts.order1_total = eta1;
    eta_parts.order2_total = eta2;
    eta_parts.order3_total = eta3;
    eta_parts.order2_increment = eta2 - eta1;
    eta_parts.order3_increment = eta3 - eta2;

    phi_parts = struct();
    phi_parts.order1_total = phi1;
    phi_parts.order2_total = phi2;
    phi_parts.order3_total = phi3;
    phi_parts.order2_increment = phi2 - phi1;
    phi_parts.order3_increment = phi3 - phi2;
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
