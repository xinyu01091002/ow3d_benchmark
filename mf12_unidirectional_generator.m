% mf12_unidirectional_generator.m
% Generate third-order unidirectional wave-group initial conditions for OceanWave3D
% using the bundled MF12 direct reconstruction interface.

clc;
clear;
close all;

CFG = struct();

% -------------------- Sea state --------------------
CFG.g = 9.81;
CFG.kp = 0.0279;
CFG.Akp_list = [0.06];
CFG.Alpha_list = [8];
CFG.kd_list = [5];
CFG.phases_deg = 0:90:270;
CFG.single_case_only = true;

% -------------------- Unidirectional spectrum --------------------
CFG.kw_left = 0.004606;
CFG.max_components = 300;
CFG.energy_keep_frac = 0.995;

% -------------------- Domain / timing --------------------
CFG.t_init_periods = -20; % Initial-condition time relative to focus in units of Tp.
CFG.t_focus = 0.0;
CFG.duration_periods = 20; % Total OW3D duration in Tp after initialization.
CFG.steps_per_period = 30;
CFG.Lx_lambda = 68; % Physical domain size in x: Lx = Lx_lambda * lambda_p.
CFG.Nx = 4097; % Number of x-grid points written to OW3D.init.
CFG.focus_x_fraction = 0.5; % Focus point as a fraction of Lx.

% -------------------- Output --------------------
CFG.output_dir = fullfile('uni initial condition', 'test_generator');
CFG.store_surface_stride = 1; % Section 8 surface output only. 1 saves every step, 2 every second step, etc.
CFG.surface_format = 1; % Keep the existing OW3D surface-output format used in this project.

setup_mf12_paths();

lambda_p = 2 * pi / CFG.kp;
Lx = CFG.Lx_lambda * lambda_p;
dx = Lx / (CFG.Nx - 1);
x = (0:CFG.Nx-1) * dx;
y = 0;

fprintf('MF12 unidirectional generator\n');
fprintf('Domain: Lx=%.3f m, Nx=%d\n', Lx, CFG.Nx);
fprintf('Grid: dx=%.3f m\n', dx);

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
    Tp = 2 * pi / omega_p;
    t_eval = CFG.t_init_periods * Tp;
    dt = Tp / CFG.steps_per_period;
    n_steps = round((CFG.duration_periods - CFG.t_init_periods) * CFG.steps_per_period);
    x_focus = CFG.focus_x_fraction * Lx;

    for Akp = Akp_list
        for Alpha = Alpha_list
            [a, b, kx, ky, meta] = build_unidirectional_spectrum(CFG, Akp, Alpha, h, x_focus);

            fprintf('kd=%.2f, Akp=%.3f, Alpha=%.1f: retained %d components\n', ...
                kd, Akp, Alpha, numel(kx));

            for phi_shift_deg = phase_list
                phase_shift = deg2rad(phi_shift_deg);
                a_shift = a * cos(phase_shift) - b * sin(phase_shift);
                b_shift = a * sin(phase_shift) + b * cos(phase_shift);

                coeffs = mf12_direct_coefficients(3, CFG.g, h, a_shift, b_shift, kx, ky, 0, 0);
                [eta, phi_surface] = mf12_direct_surface(3, coeffs, x, y, t_eval);

                eta = eta(:);
                phi_surface = phi_surface(:);

                write_path = fullfile(pwd, CFG.output_dir, ...
                    sprintf('T_init%d_Tp_Alpha_%.1f_Akp_%03d_kd%.1f_phi_%d', ...
                    CFG.t_init_periods, Alpha, round(Akp * 100), kd, phi_shift_deg));

                ensure_dir(write_path);

                write_ow3d_init(fullfile(write_path, 'OceanWave3D.init'), ...
                    eta, phi_surface, dx, Lx, dt, Akp, phase_shift);

                write_ow3d_inp(fullfile(write_path, 'OceanWave3D.inp'), ...
                    CFG, Lx, h, n_steps, dt);

                write_case_readme(fullfile(write_path, 'OW_readme.txt'), ...
                    CFG, meta, Akp, Alpha, kd, h, Tp, dx, t_eval, phi_shift_deg);

                save_visualizations(write_path, x, eta, phi_surface);

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

function [a, b, kx, ky, meta] = build_unidirectional_spectrum(CFG, Akp, Alpha, h, x_focus)
    kx_dense = linspace(CFG.kp / 400, 10 * CFG.kp, max(800, 4 * CFG.max_components));
    dk = kx_dense(2) - kx_dense(1);

    kw_right = sqrt(CFG.kp^2 / (2 * log(10^Alpha)));
    S = zeros(size(kx_dense));
    mask_low = kx_dense < CFG.kp;
    S(mask_low) = exp(-(kx_dense(mask_low) - CFG.kp).^2 / (2 * CFG.kw_left^2));
    S(~mask_low) = exp(-(kx_dense(~mask_low) - CFG.kp).^2 / (2 * kw_right^2));

    [S_sorted, idx_sort] = sort(S, 'descend');
    cum_energy = cumsum(S_sorted);
    n_keep = find(cum_energy >= CFG.energy_keep_frac * cum_energy(end), 1, 'first');
    n_keep = min(n_keep, CFG.max_components);
    idx = sort(idx_sort(1:n_keep));

    kx = kx_dense(idx);
    ky = zeros(size(kx));

    amp = S(idx) * dk;
    amp = amp * ((Akp / CFG.kp) / max(sum(amp), eps));

    omega = sqrt(CFG.g * kx .* tanh(h * kx));
    phase_focus = -kx * x_focus + omega * CFG.t_focus;
    a = amp .* cos(phase_focus);
    b = amp .* sin(phase_focus);

    meta = struct();
    meta.n_components = numel(kx);
    meta.energy_keep_frac = CFG.energy_keep_frac;
    meta.kmin = min(kx);
    meta.kmax = max(kx);
    meta.x_focus = x_focus;
end

function ensure_dir(path_str)
    if ~exist(path_str, 'dir')
        mkdir(path_str);
    end
end

function write_ow3d_init(file_name, eta, phi_surface, dx, Lx, dt, Akp, phase_shift)
    nx = numel(eta);
    ny = 1;
    dy = 10;
    Ly = ny * dy;

    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, ' H=%f nx=%d ny=%d dx=%f dy=%f akp=%f shift=%f', ...
        max(eta), nx, ny, dx, dy, Akp, phase_shift / pi);
    fprintf(f, '\n%12e %12e %d %d %12e', Lx, Ly, nx, ny, dt);

    for rx = 1:nx
        fprintf(f, '\n%12e %12e', eta(rx), phi_surface(rx));
    end
end

function write_ow3d_inp(file_name, CFG, Lx, h, n_steps, dt)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'Generated %s <-\n', datestr(now, 0));
    fprintf(f, '-1 2 <-\n');
    fprintf(f, '%d 0. %d %d 1 25 0 0 1 1 1 1 <-\n', round(Lx), h, CFG.Nx);
    fprintf(f, '3 0 3 1 1 1 <-\n');
    fprintf(f, '%d %f 1 0. 1 <-\n', round(n_steps) + 1, dt);
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

function write_case_readme(file_name, CFG, meta, Akp, Alpha, kd, h, Tp, dx, t_eval, phi_shift_deg)
    f = fopen(file_name, 'w');
    if f < 0
        error('Unable to open file for writing: %s', file_name);
    end
    cleanup = onCleanup(@() fclose(f));

    fprintf(f, 'MF12 third-order unidirectional initialization\n');
    fprintf(f, 'Generated: %s\n', datestr(now, 0));
    fprintf(f, 'Akp=%.4f, Alpha=%.1f, kp=%.5f, kd=%.2f, h=%.4f\n', Akp, Alpha, CFG.kp, kd, h);
    fprintf(f, 'Tp=%.6f s, t_init=%.6f s, phase shift=%d deg\n', Tp, t_eval, phi_shift_deg);
    fprintf(f, 'Components=%d, energy keep frac=%.5f\n', meta.n_components, meta.energy_keep_frac);
    fprintf(f, 'k-range=[%.6f, %.6f] 1/m\n', meta.kmin, meta.kmax);
    fprintf(f, 'Focus position: x=%.6f m\n', meta.x_focus);
    fprintf(f, 'Grid: Nx=%d, dx=%.6f m\n', CFG.Nx, dx);
    fprintf(f, 'Total OW3D duration after initialization: %.2f Tp\n', CFG.duration_periods);
    fprintf(f, 'Surface output stride: every %d time step(s)\n', abs(CFG.store_surface_stride));
    fprintf(f, 'Kinematic output: disabled\n');
end

function save_visualizations(write_path, x, eta, phi_surface)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1200 500]);
    tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile;
    plot(x, eta, 'b-', 'LineWidth', 1.4);
    grid on;
    xlabel('x (m)');
    ylabel('\eta (m)');
    title('Unidirectional MF12 surface elevation');

    nexttile;
    plot(x, phi_surface, 'r-', 'LineWidth', 1.4);
    grid on;
    xlabel('x (m)');
    ylabel('\phi_s (m^2/s)');
    title('Unidirectional MF12 free-surface potential');

    exportgraphics(fig, fullfile(write_path, 'mf12_unidirectional_overview.png'), 'Resolution', 200);
    close(fig);
end
