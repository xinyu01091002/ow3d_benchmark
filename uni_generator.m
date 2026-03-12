% uni_generator.m
% This code is to generate the unidirectional focused wave train initial condition, 
% with second order Dalzell and higher order VWA corrections.
% NEW: Includes Third Order Dispersion Corrections (MF12) for Loop, Potential, and Amplitude.

clc; clear; close all;

%% Configuration
CONFIG = struct();
CONFIG.L = 4097;           % Final array length
CONFIG.kp = 0.0279;        % Fixed peak wavenumber
CONFIG.g = 9.81;           % Gravity
CONFIG.Akp_list = [0.06];   % Akp values to iterate
CONFIG.Alpha_list = [1,8];  % Alpha values to iterate
CONFIG.kd_list = [0.5, 1, 5]; % Relative depth kp*h list
CONFIG.duration_periods = 20; % Total simulation duration in multiples of Tp
CONFIG.phases = 0:90:270;   % Phase shifts
CONFIG.ORDER = 5;           % Max order to compute (VWA)
CONFIG.t_init = -20;        % Initial time factor
CONFIG.output_dir = fullfile("uni initial condition", "test_generator"); % Output directory

% Constants
kp = CONFIG.kp;
lambda = 2 * pi / kp;
Nx = 68 * lambda;
dx = lambda / 30;
N_grid = Nx / dx;
kmax = 2 * pi / dx / 2;
dk = kmax / N_grid;
k_vec = dk:dk:kmax;

% Add VWA path
addpath(genpath('test functions for VWA Opensource'));
addpath(genpath('irregularWavesMF12')); % Assuming coeffsMF12 might be needed, but we inline the correction logic

% Setup Figures
figs.etah = figure('Name', 'Elevation');
figs.poth = figure('Name', 'Potential');

%% Main Loop
for target_kd = CONFIG.kd_list
    
    % --- Determine Physics for current kd ---
    h = target_kd / kp;
    wp = sqrt(CONFIG.g * kp * tanh(kp * h));
    CONFIG.Tp = 2 * pi / wp; 
    
    % Group Velocity (Linear)
    cp = sqrt(CONFIG.g / kp * tanh(kp * h));
    cw = 0.5 * cp * (1 + 2 * kp * h / sinh(2 * kp * h));
    
    fprintf('=== Physics Setup ===\n');
    fprintf('Tp: %.2fs, kd: %.2f, Depth(h): %.2fm\n', CONFIG.Tp, target_kd, h);
    fprintf('kp: %.4f, lambda: %.2fm\n', kp, lambda);
    
    for Akp = CONFIG.Akp_list
        A = Akp / kp;
        
        for Alpha = CONFIG.Alpha_list
            t = CONFIG.t_init * CONFIG.Tp;
            
            % Duration determined by user config
            T_total = (CONFIG.duration_periods - CONFIG.t_init) * CONFIG.Tp;
            dt = CONFIG.Tp / 30;
            N_steps = T_total / dt;
            
            % Vectorized Spectrum Calculation
            kw = 0.004606; 
            S = zeros(size(k_vec));
            mask_low = k_vec < kp;
            S(mask_low) = exp(-(k_vec(mask_low) - kp).^2 / (2 * kw^2));
            kw2 = sqrt(kp^2 / (2 * log(10^Alpha)));
            S(~mask_low) = exp(-(k_vec(~mask_low) - kp).^2 / (2 * kw2.^2));
            
            % Amplitude spectrum (Linear)
            S_sum = sum(S * dk);
            a = A * S * dk / S_sum;
            
            % --- MF12 Third Order Dispersion & Potential Correction ---
            % We perform this correction ONCE per sea state (A, Alpha, kd)
            fprintf('Computing MF12 3rd Order Dispersion Corrections...\n');
            [omega_corr, mu_total, a_corr] = compute_mf12_correction(a, k_vec, h, CONFIG.g);
            
            % Use corrected frequency
            w = omega_corr; 
            
            for phi_shift = CONFIG.phases
                fprintf('Processing: kd=%.1f, Akp=%.3f, Alpha=%.1f, Phi=%d\n', target_kd, Akp, Alpha, phi_shift);
                
                % 1. Linear Signal Generation (Modified with MF12 Dispersion)
                % Phase propagation using CORRECTED omega
                phi = -w * t + k_vec * t * cw + deg2rad(phi_shift);
                
                % Amplitude vector (Use corrected a_corr or original a?)
                % MF12 suggests eta = (a + F13) cos... 
                % We will use a_corr for the "Linear" input to ensure consistency
                F0 = a_corr .* numel(a_corr) .* exp(1i .* phi);
                F = complex([0, real(F0), fliplr(real(F0))], ...
                            [0, imag(F0), -fliplr(imag(F0))]);
                XX_linear = ifft(F);
                XX_linear = fftshift(XX_linear);
                
                % 2. Third Order MF12 Potential Generation
                % Instead of analytical extrapolation, we use the MF12 corrected transfer function (mu_total)
                % mu_total includes linear mu + muStar corrections
                % Note: MF12 mu is typically negative (-g/w). uni_generator expects positive potential logic.
                % phi ~ sin(theta). exp(i(theta-pi/2)) gives sin.
                % Linear potential amplitude is normally (g/w)*a.
                % MF12 mu is approx -(g/w).
                % So we use -mu_total to get positive potential amplitude coefficient.
                
                mu_eff = -mu_total; 
                
                % Potential Fourier Coeffs
                % We apply the potential transfer function to the amplitudes
                % F_pot = (a_corr * mu_eff) * exp(i(phi - pi/2))
                % Wait, should we use a_corr for potential too? Yes, consistent 3rd order.
                
                F0_pot = (a_corr .* mu_eff) .* numel(a_corr) .* exp(1i .* (phi - pi/2));
                F_p = complex([0, real(F0_pot), fliplr(real(F0_pot))], ...
                              [0, imag(F0_pot), -fliplr(imag(F0_pot))]);
                
                Phi_MF12 = ifft(F_p);
                Phi_MF12 = fftshift(Phi_MF12);
                
                % 3. Second Order (Dalzell)
                % Dalzell is a 2nd order correction. 
                % Combining MF12 (which handles 3rd order dispersion) with Dalzell (bound waves)
                % might seem redundant if we had full MF12, but here we just corrected the "Linear" mode.
                % So we still need bound waves.
                
                x_vec_unpadded = linspace(0, dx * numel(XX_linear), numel(XX_linear));
                [Dalzell_eta22, Dalzell_eta20, Dalzell_phi22, Dalzell_phi20] = ...
                    dalzell_2d(XX_linear, x_vec_unpadded, h);
                
                % Tapering (Windowing)
                N_len = length(Dalzell_eta22);
                taper_win = ones(size(Dalzell_eta22));
                taper_win(1:round(0.05*N_len)) = 0;
                taper_win(round(0.95*N_len):end) = 0;
                
                Dalzell_eta22 = Dalzell_eta22 .* taper_win;
                
                % Phi tapering
                N_phi = length(Dalzell_phi22); 
                taper_win_phi = ones(size(Dalzell_phi22));
                taper_win_phi(1:round(0.1*N_phi)) = 0;
                taper_win_phi(round(0.9*N_phi):end) = 0;
                
                XX_2_unpadded = XX_linear + Dalzell_eta22 + Dalzell_eta20;
                
                % For Potential, we use Phi_MF12 (which replaces Phi_linear) + Dalzell corrections
                PP_2_unpadded = Phi_MF12 + Dalzell_phi22  + Dalzell_phi20;
                PP_2_unpadded = PP_2_unpadded .* taper_win_phi; 
                
                % 4. Higher Order (VWA)
                vwa_x = (0:numel(XX_linear)-1)' * dx;
                
                % Call vwa_compute
                vwa_out = vwa_compute(XX_linear(:), vwa_x, h, CONFIG.g, ...
                    'nList', [3 4 5], ...
                    'analytic_side', 'neg', ...
                    'compute_eta', true, ...
                    'compute_phi_s', true);
                
                XX_final_unpadded = XX_2_unpadded;
                PP_final_unpadded = PP_2_unpadded;
                
                if CONFIG.ORDER >= 3
                    XX_final_unpadded = XX_final_unpadded + vwa_out.eta{3}.';
                    PP_final_unpadded = PP_final_unpadded + vwa_out.phi_s{3}.';
                end
                if CONFIG.ORDER >= 4
                    XX_final_unpadded = XX_final_unpadded + vwa_out.eta{4}.';
                    PP_final_unpadded = PP_final_unpadded + vwa_out.phi_s{4}.';
                end
                if CONFIG.ORDER >= 5
                    XX_final_unpadded = XX_final_unpadded + vwa_out.eta{5}.';
                    PP_final_unpadded = PP_final_unpadded + vwa_out.phi_s{5}.';
                end
                
                % 5. Padding & Output
                eta_final = asymmetric_padding(CONFIG.L, XX_final_unpadded);
                pot_final = asymmetric_padding(CONFIG.L, PP_final_unpadded);
                
                % Visualize
                plot_results(figs, dx, eta_final, pot_final, CONFIG.ORDER);
                
                % Export
                export_data(CONFIG, Akp, Alpha, CONFIG.t_init, phi_shift, eta_final, ...
                    pot_final, dx, A, kp, lambda, cw, dt, N_steps, h);
            end
        end
    end
end


%% Helper Functions for MF12 Correction

function [omega_corr, mu_total, a_corr] = compute_mf12_correction(a, k, h, g)
    % Calculates the 3rd order frequency dispersion and potential transfer function corrections
    % Based on coeffsMF12_superharmonic.m (Madsen & Fuhrman 2012)
    
    N = length(a);
    kappa = k; 
    
    % Linear quantities
    omega1 = sqrt(g*kappa.*tanh(h*kappa)); 
    F = -omega1./(kappa.*sinh(h*kappa)); 
    mu = F.*cosh(h*kappa); 
    
    % 2nd Order Terms (needed for 3rd order)
    G_2 = 1/2*h*kappa.*(2 + cosh(2*h*kappa)).*coth(h*kappa)./(sinh(h*kappa).^2); 
    F_2 = -3/4*h*omega1./(sinh(h*kappa).^4); 
    kappa_2 = 2*kappa;

    % 3rd Order Self Terms
    Upsilon = omega1.*kappa.*(-13 + 24*cosh(2*h*kappa) + cosh(4*h*kappa))./(64*sinh(h*kappa).^5); 
    Xi = 1/(4*h).*(omega1.*G_2 + F_2.*kappa_2.*sinh(h*kappa_2) - g*h*kappa.^2./(2*omega1)); 
    Omega_term = (8 + cosh(4*h*kappa))./(16*sinh(h*kappa).^4); 
    
    c = a; % Assuming a is amplitude (b=0)
    
    F13 = c.^2 .* Upsilon;
    muStar = c.^2 .* Xi;
    omega3 = c.^2 .* kappa.^2 .* Omega_term;
    
    % Interaction Terms (Double Sum Corrected)
    % Only computing interactions that affect the fundamental frequency/potential
    % (i.e. Resonant/Dispersion terms).
    % In MF12 Irregular, Eq 3.73, 3.84, 3.66 involve double summation over m.
    
    % Pre-allocate accumulated corrections
    F13_sum = zeros(size(a));
    muStar_sum = zeros(size(a));
    omega3_sum = zeros(size(a));
    
    % Speed up: Only loop over significant components? 
    % For N ~ 2000, N^2 ~ 4e6 ops, which is fine in MATLAB if inner ops are simple.
    
    % We need kx, ky. For 1D: kx = k, ky = 0.
    kx = k; ky = zeros(size(k));
    
    for n = 1:N
        % Skip if amplitude is negligible to speed up?
        if a(n) < 1e-6 * max(a), continue; end
        
        for m = [1:n-1, n+1:N]
            if a(m) < 1e-6 * max(a), continue; end
            
            % Sum interactions (n+m) --> pm=1
            pm=1;
            omega_npm = omega1(n) + pm*omega1(m);
            kx_npm = kx(n)+pm*kx(m); ky_npm = ky(n)+pm*ky(m);
            kappa_npm = abs(kx_npm); % ky is 0
            
            alpha_npm = omega_npm*cosh(h*kappa_npm);
            gamma_npm = kappa_npm*sinh(h*kappa_npm);
            beta_npm = omega_npm^2*cosh(h*kappa_npm) - g*kappa_npm*sinh(h*kappa_npm);
            
            F_npm = Gamma2_local(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                omega_npm,beta_npm, g,h);
            G_npm = Lambda2_local(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                omega_npm,alpha_npm,gamma_npm,beta_npm, g,h);
            
            % Diff interactions (n-m) --> pm=-1
            pm=-1;
            omega_nmm = omega1(n) + pm*omega1(m);
            kx_nmm = kx(n)+pm*kx(m); ky_nmm = ky(n)+pm*ky(m);
            kappa_nmm = abs(kx_nmm);
            
            alpha_nmm = omega_nmm*cosh(h*kappa_nmm);
            gamma_nmm = kappa_nmm*sinh(h*kappa_nmm);
            beta_nmm = omega_nmm^2*cosh(h*kappa_nmm) - g*kappa_nmm*sinh(h*kappa_nmm);
            
            F_nmm = Gamma2_local(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                omega_nmm,beta_nmm, g,h);
            G_nmm = Lambda2_local(omega1(n),kx(n),ky(n),kappa(n), pm*omega1(m),pm*kx(m),pm*ky(m),kappa(m), ...
                omega_nmm,alpha_nmm,gamma_nmm,beta_nmm, g,h);  
            
            % Calculate Kernels
            Ups_nm = Upsilon_nm_local(omega1(n),kx(n),ky(n),kappa(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_npm,F_nmm,G_npm,G_nmm, kappa_npm,kappa_nmm ,g,h);
            
            XI_nm = Xi_nm_local(omega1(n),kappa(n), omega1(m), G_npm,F_npm,gamma_npm, G_nmm,F_nmm,gamma_nmm, h,g);
            
            Om_nm = Omega_nm_local(omega1(n),kx(n),ky(n), omega1(m),kx(m),ky(m),kappa(m), ...
                F_npm,F_nmm,G_npm,G_nmm, kappa_npm,kappa_nmm, g,h);
            
            % Aggregate
            F13_sum(n) = F13_sum(n) + c(m)^2 * Ups_nm;
            muStar_sum(n) = muStar_sum(n) + c(m)^2 * XI_nm;
            omega3_sum(n) = omega3_sum(n) + c(m)^2 * kappa(m)^2 * Om_nm;
        end
    end
    
    % Final corrected values
    omega_corr = omega1 + (omega3 + omega3_sum) .* omega1; % Eq 3.72
    muStar = muStar + muStar_sum;
    muStar = muStar + (F13 + F13_sum) .* cosh(h*kappa); % Eq. 3.84 full
    
    mu_total = mu + muStar;
    
    F13_total = F13 + F13_sum;
    a_corr = a + F13_total; % Updated linear amplitude for free surface
    
end

% --- MF12 Local Helper Functions (Inlined) ---

function out = Lambda2_local(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,alpha_npm,gamma_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy; 
    out = h/(2*omega1n*omega1m*beta_npm)*(g*alpha_npm*(omega1n*(kappam^2 + knkm) ...
                                                     + omega1m*(kappan^2 + knkm)) ...
        + gamma_npm*(g^2*knkm + omega1n^2*omega1m^2 - omega1n*omega1m*omega_npm^2));
end

function out = Gamma2_local(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, omega_npm,beta_npm,g,h)
    knkm = knx*kmx + kny*kmy;
    out = h/(2*omega1n*omega1m*beta_npm)*(omega1n*omega1m*omega_npm*(omega_npm^2 - omega1n*omega1m) ...
        - g^2*omega1n*(kappam^2 + 2*knkm) - g^2*omega1m*(kappan^2 + 2*knkm));
end

function out = Upsilon_nm_local(omega1n,knx,kny,kappan, omega1m,kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy; 
    out = g/(4*omega1n*omega1m*cosh(h*kappan))*(omega1m*(kappan^2 - kappam^2) - omega1n*knkm) ...
        + (Gnpm + Gnmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*(g^2*knkm + omega1m^3*omega1n) ...
        - 1/(4*h*cosh(h*kappan))*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        + g*Fnpm*cosh(h*kappanpm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n + omega1m)*(knkm + kappam^2) - omega1m*kappanpm^2) ...
        + g*Fnmm*cosh(h*kappanmm)/(4*h*omega1n^2*omega1m*cosh(h*kappan))*((omega1n - omega1m)*(knkm - kappam^2) - omega1m*kappanmm^2);
end

function out = Xi_nm_local(omega1n,kappan, omega1m, Gnpm,Fnpm,gamma_npm, Gnmm,Fnmm,gamma_nmm, h,g)
    out = 1/(2*h)*(omega1m*(Gnpm - Gnmm) + Fnpm*gamma_npm + Fnmm*gamma_nmm - g*h*kappan^2/(2*omega1n));
end

function out = Omega_nm_local(omega1n,knx,kny, omega1m,kmx,kmy,kappam, Fnpm,Fnmm,Gnpm,Gnmm, kappanpm,kappanmm, g,h)
    knkm = knx*kmx + kny*kmy; 
    
    out = 1/(kappam^2)*((2*omega1m^2 + omega1n^2)/(4*omega1n*omega1m)*knkm + 1/4*kappam^2) ...
        + (Gnpm + Gnmm)/(kappam^2)*(g*knkm/(4*h*omega1n*omega1m) - omega1m^2/(4*g*h)) ...
        + omega1n/(4*g*h*kappam^2)*(Fnpm*kappanpm*sinh(h*kappanpm) + Fnmm*kappanmm*sinh(h*kappanmm)) ...
        - Fnpm*cosh(h*kappanpm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n - omega1m)*(kappam^2 + knkm) + omega1m*kappanpm^2) ...
        + Fnmm*cosh(h*kappanmm)/(4*h*omega1n*omega1m*kappam^2)*((omega1n + omega1m)*(kappam^2 - knkm) - omega1m*kappanmm^2);
end

function plot_results(figs, dx, eta, pot, order)
    xx = (1:numel(eta)) * dx;
    xx = xx - xx(end)/2;
    
    figure(figs.etah); hold on;
    plot(xx, eta, 'DisplayName', sprintf('Order %d', order));
    title('Surface Elevation'); grid on; box on;
    
    figure(figs.poth); hold on;
    plot(xx, pot, 'DisplayName', sprintf('Order %d', order));
    title('Velocity Potential'); grid on; box on;
    drawnow;
end

function export_data(CFG, Akp, Alpha, t_init, phi_shift, eta_S22, pot_S22, dx, A, kp, lambda, cw, dt, N_steps, h)
    % Prepare directory
    write_path = fullfile(pwd, CFG.output_dir, sprintf("T_init%d_Tp_Alpha_%.1f_Akp_%.3d_kd%.1f_phi_%d", ...
        t_init, Alpha, round(Akp*100), kp*h, phi_shift));
        
    if ~exist(write_path, 'dir')
        mkdir(write_path);
        disp(['Directory created: ', write_path]);
    else
        % disp('Directory exists.');
    end
    
    % 1. OceanWave3D.init
    file_name = fullfile(write_path, 'OceanWave3D.init');
    nx = length(eta_S22);
    ny = 1;
    dy = 10; % arbitrary y-width
    Lx = nx * dx; % Or based on max-min
    Ly = ny * dy;
    
    f = fopen(file_name, 'w');
    fprintf(f, ' H=%f nx=%d ny=%d dx=%f dy=%f akp=%f shift=%f', ...
        max(eta_S22), nx, ny, dx, dy, Akp, phi_shift/pi);
    fprintf(f, '\n%12e %12e %d %d %12e', Lx, Ly, nx, ny, dt);
    
    % Write data (row-major implicit in loop?)
    for ry = 1:ny
        for rx = 1:nx
            fprintf(f, '\n%12e %12e', eta_S22(rx), pot_S22(rx));
        end
    end
    fclose(f);
    
    % 2. OceanWave3D.inp
    file_name = fullfile(write_path, 'OceanWave3D.inp');
    f = fopen(file_name, 'w');
    fprintf(f, 'Generated %s <-\n', datestr(now, 0));
    fprintf(f, '-1 2 <-\n');
    fprintf(f, '%d 0. %d %d 1 25 0 0 1 1 1 1 <-\n', round(CFG.L*dx), h, CFG.L);
    fprintf(f, '3 0 3 1 1 1 <-\n');
    fprintf(f, '%d %f 1 0. 1 <-\n', round(N_steps)+1, dt);
    fprintf(f, '9.81 <-\n');
    fprintf(f, '1 3 0 55 1e-6 1e-6 1 V 1 1 20 <-\n');
    fprintf(f, '0.05 1.00 1.84 2 0 0 1 6 32 <-\n');
    fprintf(f, '8 1 <-\n');
    fprintf(f, '1 0 <-\n');
    fprintf(f, '0 6 10 0.08 0.08 0.4 <-\n');
    fprintf(f, '0 8. 3 X 0.0 <-\n');
    fprintf(f, '0 0 <-\n');
    fprintf(f, '0 2.0 2 0 0 1 0 <-\n');
    fprintf(f, '0 <-\n');
    fprintf(f, '33  8. 2. 80. 20. -1 -11 100. 50. run06.el 22.5 1.0 3.3 <-\n');
    fclose(f);
    
    % 3. Readme
    file_name = fullfile(write_path, 'OW_readme.txt');
    f = fopen(file_name, 'w');
    t_val = t_init * CFG.Tp; 
    fprintf(f, 'Data for Tim T=%1.2fTp H=%f dx=%f, dt=%f\n', t_val/CFG.Tp, max(eta_S22), dx, dt);
    fprintf(f, 'A=%f, Akp=%f, kp=%f, lambda=%f, Tp=%f\n', A, Akp, kp, lambda, CFG.Tp);
    fprintf(f, 'lambda=%f dx, T=%f dt, cw=%f, CFL=%f, kd=%f\n', ...
        lambda/dx, CFG.Tp/dt, cw, cw*dt/dx, kp*h);
    fprintf(f, 'Last modified in %s', datestr(now, 0));
    fclose(f);
end

function padded_array = asymmetric_padding(L, array)
    % This code is to perform asymmertic zero padding to the original array
    totalPadding = L - numel(array);
    if totalPadding < 0
        warning('Output length L is smaller than input array.');
        padded_array = array;
        return;
    end
    
    leftRatio = 1;
    rightRatio = 20;
    totalRatio = leftRatio + rightRatio;
    
    leftPad = round(totalPadding * (leftRatio / totalRatio));
    rightPad = totalPadding - leftPad;
    
    padded_array = padarray(array, [0, leftPad], 0, 'pre');
    padded_array = padarray(padded_array, [0, rightPad], 0, 'post');
    
    % Smooth transition
    % Left side
    for i = leftPad:-1:2
        padded_array(i-1) = padded_array(i) / 2;
    end
    
    % Right side
    start_idx = leftPad + length(array);
    for i = 1:rightPad
        padded_array(start_idx + i) = padded_array(start_idx + i - 1) / 2;
    end
end
