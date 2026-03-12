%% Test Script for Universal VWA Calculation Function (vwa_compute.m)
% Comparison with OW3D (Exact/FNPF) Reference Data
% This script loads an existing dataset from processed_eta22 and compares
% the VWA bound harmonic prediction against the ground truth.

clc; clear; close all;

%% 1. Path and Configuration
% Locate the processed_eta22 folder (sibling to current folder)
current_dir = fileparts(mfilename('fullpath'));
% Assuming structure: .../Directional/test functions/test_script.m
% We need to go up to Directional/processed_eta22
project_dir = fileparts(current_dir); 
data_dir = fullfile(project_dir, 'processed_eta33');
export_folder = fullfile(current_dir, 'figures_comparison');
if ~isfolder(export_folder), mkdir(export_folder); end

% Configuration (matching available data file)
config = struct();
config.dataset_label = 'test6';
config.kd = 5.0;
config.spread_deg = 15;
config.Akp = 0.12;
config.time_step = 800;

% Construct filename
snapshot_tag = sprintf('OW3D_eta33_%s_kd%.1f_spread%d_Akp%.2f_t%05d.mat', ...
    config.dataset_label, config.kd, config.spread_deg, config.Akp, config.time_step);
snapshot_path = fullfile(data_dir, snapshot_tag);


fprintf('Loading data from: %s\n', snapshot_path);
S = load(snapshot_path);

%% 2. Extract Data
fprintf('Loading eta22 data from: %s\n', snapshot_path);
S2 = load(snapshot_path);
eta11_exact = S2.eta11;  % [Nx x Ny]
eta22_exact = S2.eta22;  % [Nx x Ny]
X = S2.X;
Y = S2.Y;

x_vec = X(:,1);
y_vec = Y(1,:);
Nx = length(x_vec);
Ny = length(y_vec);
dx = x_vec(2) - x_vec(1);

% Try to load eta33/eta44/eta55 if available
snapshot_tag33 = strrep(snapshot_tag, 'eta22', 'eta33');
snapshot_path33 = fullfile(data_dir, '../processed_eta33', snapshot_tag33);
eta33_exact = [];
eta44_exact = [];
eta55_exact = [];

if isfile(snapshot_path33)
    fprintf('Loading higher order data from: %s\n', snapshot_tag33);
    S3 = load(snapshot_path33);
    
    if isfield(S3, 'eta33')
        eta33_exact = S3.eta33;
    elseif isfield(S3, 'eta_harmonics') && size(S3.eta_harmonics, 1) >= 3
         eta33_exact = squeeze(S3.eta_harmonics(3,:,:));
    end

    if isfield(S3, 'eta44')
        eta44_exact = S3.eta44;
    elseif isfield(S3, 'eta_harmonics') && size(S3.eta_harmonics, 1) >= 4
        % eta44 is usually index 4. This bin contains 0, 4, 8... harmonics.
        % We should filter it to isolate n=4.
        fprintf('  Attempting to reconstruct eta44 from harmonic 4...\n');
        try
             raw_harm4 = squeeze(S3.eta_harmonics(4,:,:));
             X3 = S3.X; 
             Y3 = S3.Y;
             x3_vec = X3(:,1);
             y3_vec = Y3(1,:);
             kp_val = 2 * pi / lambda;
             
             addpath(project_dir);
             if exist('frequency_filtering_2d', 'file')
                 eta44_exact = frequency_filtering_2d(raw_harm4, x3_vec, y3_vec, kp_val, 4);
                 fprintf('  Success: Extracted eta44 via frequency_filtering_2d.\n');
             else
                 eta44_exact = raw_harm4;
                 fprintf('  Warning: frequency_filtering_2d not found. Using raw harmonic 4 (may include DC).\n');
             end
        catch ME
            fprintf('  Failed to extract eta44: %s\n', ME.message);
        end
    end
    
    if isfield(S3, 'eta55')
        eta55_exact = S3.eta55;
    elseif isfield(S3, 'eta_harmonics') && size(S3.eta_harmonics, 1) >= 1
         % Extract eta55 from eta_harmonics(1) (aliased with eta11) using filtering
         % This logic mimics extract_eta33_from_OW3D.m
         fprintf('  Attempting to reconstruct eta55 from harmonic 1 (aliased)...\n');
         try
             raw_harm1 = squeeze(S3.eta_harmonics(1,:,:));
             X3 = S3.X; 
             Y3 = S3.Y;
             x3_vec = X3(:,1);
             y3_vec = Y3(1,:);
             kp_val = 2 * pi / lambda;
             
             % Ensure access to frequency_filtering_2d in parent folder
             addpath(project_dir);
             
             if exist('frequency_filtering_2d', 'file')
                 eta55_exact = frequency_filtering_2d(raw_harm1, x3_vec, y3_vec, kp_val, 5);
                 fprintf('  Success: Extracted eta55 via frequency_filtering_2d.\n');
             else
                 fprintf('  Warning: frequency_filtering_2d.m not found. Cannot extract eta55.\n');
             end
         catch ME
             fprintf('  Failed to extract eta55: %s\n', ME.message);
         end
    end
else
    fprintf('Warning: Higher order data not found at %s. Skipping comparison.\n', snapshot_path33);
end

% Physical Parameters
lambda = 225;                   % Reference wavelength
kp_ref = 2 * pi / lambda;       % ~0.0279
h = config.kd / kp_ref;         
g = 9.81;
omega_p = sqrt(g * kp_ref * tanh(kp_ref * h));

fprintf('Grid Size: %d x %d\n', Nx, Ny);
fprintf('Parameters: h = %.2f m, kp = %.4f rad/m, g = %.2f, omega_p = %.4f\n', h, kp_ref, g, omega_p);

%% 3. Configure and Run VWA
% Compute orders 2 through 5, enabling phi_s calculation
vwa_opts = { ...
    'nList', [2 3 4 5], ...           
    'analytic_side', 'pos', ... 
    'phi_take', 'imag', ...
    'compute_eta', true, ...
    'compute_phi_s', true      
};

fprintf('Running vwa_compute (Orders 2-5, eta & phi)...\n');
tic;
try
    out = vwa_compute(eta11_exact, x_vec, h, g, vwa_opts{:});
catch ME
    fprintf('Error: %s\n', ME.message);
    rethrow(ME);
end
toc;

%% 4. Comparison and Visualization
% Normalization factor for plots: use (A^n * k^(n-1)) or similar?
% Usually papers plot eta^(n) / (A^n * k^(n-1)) or simply eta^(n) / A^n * factor
% Here we use a coherent normalization: eta / A. Or just raw values?
% The previous script used eta22 / (A^2 k). 
% Let's generalize: Scaling ~ A^n.
% For n=2: A^2 k
% For n=3: A^3 k^2
% For n=4: A^4 k^3 ...
A_ref = config.Akp / kp_ref;

% Locate Centerline
mid_idx = floor(Ny/2) + 1;
x_norm = x_vec / lambda;

% Determine plotting range based on eta11 peak
eta11_center = eta11_exact(:, mid_idx);
[~, peak_idx] = max(abs(eta11_center));
peak_x_norm = x_norm(peak_idx);
plot_range_lambda = 3; 
plot_xlim = [peak_x_norm - plot_range_lambda, peak_x_norm + plot_range_lambda];

% Prepare Figure
fig_comp = figure('Name', 'VWA High Order Comparison', 'Color', 'w', 'Position', [50, 50, 1200, 800]);

orders_to_plot = [2 3 4 5];
for i = 1:length(orders_to_plot)
    n = orders_to_plot(i);
    
    % Retrieve VWA result
    if n <= length(out.eta) && ~isempty(out.eta{n})
        eta_vwa_raw = out.eta{n};
        cta_vwa_center = eta_vwa_raw(:, mid_idx);
    else
        cta_vwa_center = zeros(Nx, 1);
    end
    
    % Retrieve Exact result
    eta_exact_center = [];
    if n == 2
        eta_exact_center = eta22_exact(:, mid_idx);
    elseif n == 3 && ~isempty(eta33_exact)
        eta_exact_center = eta33_exact(:, mid_idx);
    elseif n == 4 && ~isempty(eta44_exact)
        eta_exact_center = eta44_exact(:, mid_idx);
        eta_exact_center=frequency_filtering(eta_exact_center,x_vec,0.0279,4);
    elseif n == 5 && ~isempty(eta55_exact)
        eta_exact_center = eta55_exact(:, mid_idx);
    end
    
    % Normalization: eta^(n) / (A^n * kp^(n-1))
    scale_factor = (A_ref^n) * (kp_ref^(n-1));
    % Avoid division by zero
    if scale_factor == 0, scale_factor = 1; end
    
    val_vwa = cta_vwa_center / scale_factor;
    val_exact = [];
    if ~isempty(eta_exact_center)
        val_exact = eta_exact_center / scale_factor;
        
        % Sign check (only if we have exact data to compare)
        % For higher orders, sign conventions can vary. Flip VWA if anti-correlated.
        tmp_corr = corrcoef(val_exact, val_vwa);
        if tmp_corr(1,2) < -0.5
             val_vwa = -val_vwa;
        end
    end
    
    % Plotting in Subplot
    subplot(2, 2, i);
    hold on;
    if ~isempty(val_exact)
        plot(x_norm, val_exact, 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Exact (OW3D)');
    end
    plot(x_norm, val_vwa, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'VWA');
    
    grid on; box on;
    xlim(plot_xlim);
    xlabel('$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(sprintf('$\\eta^{(%d)}$ (norm)', n), 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Order %d', n), 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    
    % Only add legend to the first plot to save space, or all if preferred
    if i == 1
        legend('Location', 'best', 'FontSize', 10, 'Interpreter', 'latex');
    end
    
    % Display max amplitudes
    max_vwa = max(abs(val_vwa));
    msg = sprintf('Max VWA: %.3f', max_vwa);
    if ~isempty(val_exact)
        max_exact = max(abs(val_exact));
        msg = [msg, sprintf(', Exact: %.3f', max_exact)];
        rmse = sqrt(mean((val_exact - val_vwa).^2));
        msg = [msg, sprintf(', RMSE: %.3f', rmse)];
    end
    text(0.05, 0.95, msg, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 2);
end

sgtitle(sprintf('Comparison of Bound Harmonics (Orders 2-5)\nDataset: %s, Spread: %d deg', config.dataset_label, config.spread_deg));

% Save Figure
saveas(fig_comp, fullfile(export_folder, 'comparison_vwa_orders_2to5.png'));
fprintf('Comparison figure saved.\n');

%% 4b. Comparison and Visualization - PHI
% Load Phi Exact Data
% Basic Loading Logic for Phi
phi22_exact = []; phi33_exact = []; phi44_exact = []; phi55_exact = [];

if isfield(S2, 'phi22'), phi22_exact = S2.phi22; end

% Phi 3,4,5 from S3
if isfile(snapshot_path33)
    if isfield(S3, 'phi33'), phi33_exact = S3.phi33; end
    if isfield(S3, 'phi_harmonics')
        if size(S3.phi_harmonics, 1) >= 4
             % Filter phi44 from harmonic 4
             raw_phi4 = squeeze(S3.phi_harmonics(4,:,:));
             addpath(project_dir);
             if exist('frequency_filtering_2d', 'file')
                 phi44_exact = frequency_filtering_2d(raw_phi4, S3.X(:,1), S3.Y(1,:), kp_ref, 4);
             else
                 phi44_exact = raw_phi4; 
             end
        end
        if size(S3.phi_harmonics, 1) >= 1
             % Filter phi55 from harmonic 1
             raw_phi1 = squeeze(S3.phi_harmonics(1,:,:));
             if exist('frequency_filtering_2d', 'file')
                 phi55_exact = frequency_filtering_2d(raw_phi1, S3.X(:,1), S3.Y(1,:), kp_ref, 5);
             end
        end
    end
end

% Phi Normalization Factor: A^n * (g / omega_p) * k^(n-1) approx.
fig_comp_phi = figure('Name', 'VWA High Order Comparison (Phi)', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

for i = 1:length(orders_to_plot)
    n = orders_to_plot(i);
    
    % VWA Phi
    if n <= length(out.phi_s) && ~isempty(out.phi_s{n})
        phi_vwa_raw = out.phi_s{n};
        phi_vwa_center = phi_vwa_raw(:, mid_idx);
    else
        phi_vwa_center = zeros(Nx, 1);
    end
    
    % Exact Phi
    phi_exact_center = [];
    if n == 2 && ~isempty(phi22_exact), phi_exact_center = phi22_exact(:, mid_idx); end
    if n == 3 && ~isempty(phi33_exact), phi_exact_center = phi33_exact(:, mid_idx); end
    if n == 4 && ~isempty(phi44_exact), phi_exact_center = phi44_exact(:, mid_idx); end
    if n == 5 && ~isempty(phi55_exact), phi_exact_center = phi55_exact(:, mid_idx); end
    
    % Normalization
    scale_factor_phi = (A_ref^n) * (g / omega_p) * (kp_ref^(n-1));
    if scale_factor_phi == 0, scale_factor_phi = 1; end
    
    val_vwa = phi_vwa_center / scale_factor_phi;
    val_exact = [];
    if ~isempty(phi_exact_center)
        val_exact = phi_exact_center / scale_factor_phi;
        % Sign check
        tmp_corr = corrcoef(val_exact, val_vwa);
        if tmp_corr(1,2) < -0.5
             val_vwa = -val_vwa;
        end
    end
    
    % Plot
    subplot(2, 2, i);
    hold on;
    if ~isempty(val_exact)
        plot(x_norm, val_exact, 'Color', [0 0.4470 0.7410], 'LineWidth', 2, 'DisplayName', 'Exact (OW3D)');
    end
    plot(x_norm, val_vwa, 'Color', [0.8500 0.3250 0.0980], 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'VWA');
    
    grid on; box on;
    xlim(plot_xlim);
    xlabel('$x / \lambda$', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(sprintf('$\\phi_s^{(%d)}$ (norm)', n), 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('Order %d ($\\phi_s$)', n), 'Interpreter', 'latex', 'FontSize', 13);
    set(gca, 'FontSize', 11, 'LineWidth', 1.2);
    
    if i == 1, legend('Location', 'best', 'FontSize', 10, 'Interpreter', 'latex'); end
    
    % Stats
    max_vwa = max(abs(val_vwa));
    msg = sprintf('Max VWA: %.3f', max_vwa);
    if ~isempty(val_exact)
        max_exact = max(abs(val_exact));
        msg = [msg, sprintf(', Exact: %.3f', max_exact)];
    end
    text(0.05, 0.95, msg, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 2);
end

sgtitle(sprintf('Comparison of Bound Potentials (Orders 2-5)\nDataset: %s, Spread: %d deg', config.dataset_label, config.spread_deg));
saveas(fig_comp_phi, fullfile(export_folder, 'comparison_vwa_phi_orders_2to5.png'));
fprintf('Phi comparison figure saved.\n');


