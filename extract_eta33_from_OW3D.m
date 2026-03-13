% extract_eta33_from_OW3D.m
% -------------------------------------------------------------------------
% Utility script to reconstruct the OW3D-bound third-order \eta_{33} field from
% the directional datasets. The workflow mirrors extract_eta22_from_OW3D and
% reuses the four-phase Hilbert reconstruction to recover the third-harmonic
% amplitude for fast comparisons.
% -------------------------------------------------------------------------
clc; clear; close all;

%% --- User configuration -------------------------------------------------
dataset_label   = 'test1';
kd_target       = 5;
spread_angle    = 5;
Akp             = 0.02;
lambda          = 225;
phi_shifts      = 0:90:270;
time_step       = 1200;
result_section  = 256;

%% --- Derived paths ------------------------------------------------------
this_dir     = fileparts(mfilename('fullpath'));
data_root    = fullfile(this_dir, dataset_label);
folder_stamp = compose('kd%.1f_spread_%d_Akp_%.2f_phi_shift_%d', ...
    kd_target, spread_angle, Akp, phi_shifts);

%% --- Read OW3D phases ---------------------------------------------------
eta_phases = [];
phi_phases = [];
for idx = 1:numel(phi_shifts)
    phase_folder = fullfile(data_root, folder_stamp{idx});
    if ~isfolder(phase_folder)
        error('Missing phase folder: %s', phase_folder);
    end
    bin_name = sprintf('EP_%05d.bin', time_step);
    bin_path = fullfile(phase_folder, bin_name);
    if ~isfile(bin_path)
        error('Missing OW3D output: %s', bin_path);
    end

    [X, Y, E, P] = ReadBinFile(bin_path); %#ok<NASGU>
    eta_phases(idx,:,:) = E;
    phi_phases(idx,:,:) = P;
    fprintf('Loaded %s', bin_path);
end

%% --- Harmonic reconstruction (eta33) ----------------------------------
four_phase_coef = [
    0.25  0  -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0  -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];
[eta_harmonics, phi_harmonics] = reconstruct_directional_harmonics( ...
    eta_phases, phi_phases, four_phase_coef);

eta11 = squeeze(eta_harmonics(1,:,:));
eta22 = squeeze(eta_harmonics(2,:,:));
eta33 = squeeze(eta_harmonics(3,:,:));
eta44 = squeeze(eta_harmonics(4,:,:));
phi22 = squeeze(phi_harmonics(2,:,:));
phi33 = squeeze(phi_harmonics(3,:,:));
phi44 = squeeze(phi_harmonics(4,:,:));
% Extract 5th harmonic from the 1st harmonic group (aliased)
kp = 2 * pi / lambda;
x_vec_full = X(:,1);
y_vec_full = Y(1,:);
eta11 = frequency_filtering_2d(eta11, x_vec_full, y_vec_full, kp, 1);
eta55 = frequency_filtering_2d(squeeze(eta_harmonics(1,:,:)), x_vec_full, y_vec_full, kp, 5);
phi55 = frequency_filtering_2d(squeeze(phi_harmonics(1,:,:)), x_vec_full, y_vec_full, kp, 5);

fprintf('eta33 amplitude range: [%.3f, %.3f] m\n', min(eta33(:)), max(eta33(:)));

%% --- Save results -------------------------------------------------------
output_dir = fullfile(this_dir, 'processed_eta33');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

meta = struct('dataset', dataset_label, 'kd', kd_target, ...
    'spread', spread_angle, 'Akp', Akp, 'time_step', time_step, ...
    'phi_shifts', phi_shifts);

save(fullfile(output_dir, sprintf( ...
    'OW3D_eta33_%s_kd%.1f_spread%d_Akp%.2f_t%05d.mat', ...
    dataset_label, kd_target, spread_angle, Akp, time_step)), ...
    'X', 'Y', 'eta11', 'eta22', 'eta33', 'eta44', 'eta55', ...
    'phi22', 'phi33', 'phi44', 'phi55', ...
    'eta_harmonics', 'phi_harmonics', 'meta');

disp('eta33 extraction complete. Results stored under processed_eta33/.');

%% --- Quick visualization ----------------------------------------------
x_vec = X(:,1) - X(end,1)/2;
y_vec = Y(1,:) - Y(1,end)/2;

figure('Color', 'w', 'Position', [120 100 920 620]);
subplot(1,2,1);
imagesc(x_vec, y_vec, eta33');
axis image; colorbar;
colormap(mymap('RdBu'));
title(sprintf('OW3D | eta_{33} | kd=%.1f, spread=%d^\circ', kd_target, spread_angle));
xlabel('x (m)'); ylabel('y (m)');
set(gca, 'FontSize', 12, 'LineWidth', 1.1, 'YDir', 'normal');

y_idx = min(max(1, result_section), size(eta33,2));
subplot(1,2,2);
plot(x_vec, eta11(:, y_idx), 'k-', 'LineWidth', 1.2, 'DisplayName', '\eta_{11}'); hold on;
plot(x_vec, eta22(:, y_idx), 'r-', 'LineWidth', 1.2, 'DisplayName', '\eta_{22}');
plot(x_vec, eta33(:, y_idx), 'b--', 'LineWidth', 1.2, 'DisplayName', '\eta_{33}');
plot(x_vec, phi33(:, y_idx), 'm-.', 'LineWidth', 1.1, 'DisplayName', '\phi_{33}');
hold off;
legend('Location', 'best'); grid on; box on;
xlabel('x (m)'); ylabel(sprintf('Section y_{%d}', y_idx));
title('Line-out cross section');
set(gca, 'FontSize', 12, 'LineWidth', 1.1);

%% --- Local helper functions --------------------------------------------
function [eta_harmonics, phi_harmonics] = reconstruct_directional_harmonics( ...
        eta_phases, phi_phases, coef)
    eta_hilbert = hilbert2d_all(eta_phases);
    phi_hilbert = hilbert2d_all(phi_phases);

    All_eta = cat(1, real(eta_phases), -imag(eta_hilbert));
    All_phi = cat(1, real(phi_phases), -imag(phi_hilbert));

    eta_harmonics = zeros(4, size(eta_phases,2), size(eta_phases,3));
    phi_harmonics = zeros(4, size(phi_phases,2), size(phi_phases,3));
    for idx = 1:4
        eta_harmonics(idx,:,:) = sum(All_eta .* reshape(coef(idx,:), [8,1,1]), 1);
        phi_harmonics(idx,:,:) = sum(All_phi .* reshape(coef(idx,:), [8,1,1]), 1);
    end
end

function X_hilbert = hilbert2d_all(X)
    [P, M, N] = size(X);
    X_hilbert = zeros(size(X));
    for p = 1:P
        FX = fft2(squeeze(X(p,:,:)));
        H = zeros(M, N);
        H(1:floor(M/2), :) = 1;
        H(ceil(M/2):end, :) = -1;
        X_hilbert(p,:,:) = ifft2(FX .* H);
    end
end

function field_out = frequency_filtering(field_in, x_vec, kp, n)
    % frequency_filtering  Filter 2D field in x-direction to keep harmonic n
    % field_in: (Nx, Ny)
    % x_vec: (Nx, 1)
    
    dx = x_vec(2) - x_vec(1);
    Nx = size(field_in, 1);
    dk = 2*pi / (Nx * dx);
    k = [0:ceil(Nx/2)-1, -floor(Nx/2):-1]' * dk;
    
    % Bandwidth
    k_target = n * kp;
    k_width = 0.4 * kp; % +/- 0.4 kp to be safe
    
    mask = (abs(abs(k) - k_target) < k_width);
    
    % Apply mask in Fourier domain (x-direction)
    field_fft = fft(field_in, [], 1);
    field_filtered_fft = field_fft .* mask;
    field_out = ifft(field_filtered_fft, [], 1);
    
    if isreal(field_in)
        field_out = real(field_out);
    end
end

