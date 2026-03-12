% mf12_directional_generator.m
% Generate third-order directional wave-group initial conditions for OceanWave3D
% using the MF12 spectral reconstruction workflow.

clc;
clear;
close all;

CFG = struct();

% -------------------- Sea state --------------------
CFG.g = 9.81;
CFG.kp = 0.0279;
CFG.Akp_list = [0.02];
CFG.Alpha_list = [8];
CFG.kd_list = [1];
CFG.phases_deg = 0:90:270;
CFG.single_case_only = true;

% -------------------- Directional wave-group definition --------------------
CFG.heading_deg = 0;
CFG.spread_deg = 5;
CFG.energy_keep_frac = 0.995;
CFG.max_components = 200;

% -------------------- Domain / timing --------------------
CFG.Tp = 12;
CFG.t_init_periods = -15;
CFG.t_focus = 0.0;
CFG.duration_periods = 30;
CFG.steps_per_period = 30;
CFG.Lx_lambda = 30;
CFG.Ly_lambda = 20;
CFG.Nx = 1025;
CFG.Ny = 257;

% -------------------- Output --------------------
CFG.output_dir = 'mf12_directional_generator';

setup_mf12_paths();

lambda_p = 2 * pi / CFG.kp;
Lx = CFG.Lx_lambda * lambda_p;
Ly = CFG.Ly_lambda * lambda_p;
dx = Lx / CFG.Nx;
dy = Ly / CFG.Ny;
dt = CFG.Tp / CFG.steps_per_period;
N_steps = round(CFG.duration_periods * CFG.steps_per_period);

fprintf('MF12 directional generator\n');
fprintf('Domain: Lx=%.3f m, Ly=%.3f m, Nx=%d, Ny=%d\n', Lx, Ly, CFG.Nx, CFG.Ny);
fprintf('Grid: dx=%.3f m, dy=%.3f m\n', dx, dy);

kd_list = CFG.kd_list;
Akp_list = CFG.Akp_list;
Alpha_list = CFG.Alpha_list;
phase_list = CFG.phases_deg;

if CFG.single_case_only
    kd_list = kd_list(1);
    Akp_list = Akp_list(1);
    Alpha_list = Alpha_list(1);
    phase_list = phase_list(1);
end

for kd = kd_list
    h = kd / CFG.kp;
    omega_p = sqrt(CFG.g * CFG.kp * tanh(CFG.kp * h));
    Tp_kd = 2 * pi / omega_p;
    t_eval = CFG.t_init_periods * Tp_kd;
    dt_kd = Tp_kd / CFG.steps_per_period;
    N_steps_kd = round(CFG.duration_periods * CFG.steps_per_period);

    for Akp = Akp_list
        for Alpha = Alpha_list
            [kx, ky, amp_base, meta] = build_directional_group_spectrum(CFG, Lx, Ly, h, Akp, Alpha);

            fprintf('kd=%.2f, Akp=%.3f, Alpha=%.1f: retained %d components\n', ...
                kd, Akp, Alpha, numel(kx));

            xf = 0.5 * Lx;
            yf = 0.5 * Ly;
            omega_lin = sqrt(CFG.g * hypot(kx, ky) .* tanh(h * hypot(kx, ky)));

            for phi_shift_deg = phase_list
                phase_focus = -(kx * xf + ky * yf) + omega_lin * CFG.t_focus + deg2rad(phi_shift_deg);
                a = amp_base .* cos(phase_focus);
                b = amp_base .* sin(phase_focus);

                coeffs_lin = mf12_spectral_coefficients(1, CFG.g, h, a, b, kx, ky, 0, 0);
                coeffs = mf12_spectral_coefficients(3, CFG.g, h, a, b, kx, ky, 0, 0);
                [eta_lin, phi_lin] = mf12_spectral_surface(coeffs_lin, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);
                [eta_mf12, phi_mf12] = mf12_spectral_surface(coeffs, Lx, Ly, CFG.Nx, CFG.Ny, t_eval);

                % mf12_spectral_surface returns arrays of size [Ny x Nx].
                % Transpose to match the historical OW3D export convention used in this repo: [Nx x Ny].
                eta_out = eta_mf12.';
                phi_out = phi_mf12.';

                write_path = fullfile(pwd, CFG.output_dir, ...
                    sprintf('T_init%d_Tp_kd%.1f_spread_%d_heading_%d_Akp_%03d_alpha_%.1f_phi_%d', ...
                    CFG.t_init_periods, kd, round(CFG.spread_deg), round(CFG.heading_deg), ...
                    round(Akp * 100), Alpha, phi_shift_deg));

                ensure_dir(write_path);

                write_ow3d_init(fullfile(write_path, 'OceanWave3D.init'), ...
                    eta_out, phi_out, dx, dy, Lx, Ly, dt_kd, Akp, phi_shift_deg);

                write_ow3d_inp(fullfile(write_path, 'OceanWave3D.inp'), ...
                    Lx, Ly, h, CFG.Nx, CFG.Ny, N_steps_kd, dt_kd);

                write_readme(fullfile(write_path, 'OW_readme.txt'), ...
                    CFG, meta, Akp, Alpha, kd, h, Tp_kd, dx, dy, t_eval, phi_shift_deg);

                save_visualizations(write_path, eta_lin, eta_mf12, phi_lin, phi_mf12, Lx, Ly);

                fprintf('  wrote: %s\n', write_path);
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

function write_ow3d_inp(file_name, Lx, Ly, h, nx, ny, n_steps, dt)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'Data for MF12 directional wave-group initialization %s <-\n', datestr(now, 0));
    fprintf(f, '-1 2 <-\n');
    fprintf(f, '%d %d %d %d %d %d 0 0 1 1 1 1 <-\n', round(Lx), round(Ly), round(h), nx, ny, 9);
    fprintf(f, '4 4 4 1 1 1 <-\n');
    fprintf(f, '%d %f 1 0. 1 <-\n', n_steps + 1, dt);
    fprintf(f, '9.81 <-\n');
    fprintf(f, '1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n');
    fprintf(f, '0.05 1.00 1.84 2 0 0 1 6 32 <-\n');
    fprintf(f, '10 1 <-\n');
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
    fprintf(f, 'Grid: Nx=%d, Ny=%d, dx=%.6f m, dy=%.6f m\n', CFG.Nx, CFG.Ny, dx, dy);
    fprintf(f, 'Model: MF12 spectral coefficients/order-3, single directional wave group\n');
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
