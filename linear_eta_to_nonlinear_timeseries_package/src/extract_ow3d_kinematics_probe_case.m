function ow3d = extract_ow3d_kinematics_probe_case(case_cfg)
%EXTRACT_OW3D_KINEMATICS_PROBE_CASE Extract eta/phi/u/w probe time series from four kinematics files.

required = {'data_root', 'folder_pattern', 'kinematics_file_id', 'probe_sigma_mode'};
for k = 1:numel(required)
    if ~isfield(case_cfg, required{k})
        error('case_cfg must contain field ''%s''.', required{k});
    end
end
if ~isfield(case_cfg, 'probe_x_mode') || isempty(case_cfg.probe_x_mode)
    case_cfg.probe_x_mode = 'fixed_index';
end

phi_shifts = 0:90:270;
data_by_phase = cell(1, numel(phi_shifts));

for idx = 1:numel(phi_shifts)
    case_folder = fullfile(case_cfg.data_root, sprintf(case_cfg.folder_pattern, phi_shifts(idx)));
    if case_cfg.kinematics_file_id < 10
        file_name = sprintf('Kinematics0%d.bin', case_cfg.kinematics_file_id);
    else
        file_name = sprintf('Kinematics%d.bin', case_cfg.kinematics_file_id);
    end
    kin_path = fullfile(case_folder, file_name);
    if ~isfile(kin_path)
        error('Missing kinematics file: %s', kin_path);
    end
    data_by_phase{idx} = read_ow3d_kinematics_file(kin_path, 'uncorrected');
end

ref = data_by_phase{1};
sigma_vec = ref.sigma(:);
switch lower(case_cfg.probe_sigma_mode)
    case 'surface'
        [~, sigma_idx] = max(sigma_vec);
    case 'index'
        sigma_idx = case_cfg.probe_sigma_index;
    otherwise
        error('Unsupported probe_sigma_mode: %s', case_cfg.probe_sigma_mode);
end

Nt = numel(ref.t);
Nx = size(ref.eta, 2);
raw_eta_all = zeros(Nt, 4, Nx);
raw_phi_all = zeros(Nt, 4, Nx);
raw_u_all = zeros(Nt, 4, Nx);
raw_w_all = zeros(Nt, 4, Nx);

for idx = 1:numel(phi_shifts)
    d = data_by_phase{idx};
    raw_eta_all(:, idx, :) = squeeze(d.eta(:, :, 1));
    raw_phi_all(:, idx, :) = squeeze(d.phi(:, sigma_idx, :, 1));
    raw_u_all(:, idx, :) = squeeze(d.u(:, sigma_idx, :, 1));
    raw_w_all(:, idx, :) = squeeze(d.w(:, sigma_idx, :, 1));
end

[x_idx, probe_info] = resolve_probe_x_index_local(case_cfg, raw_eta_all, ref.x(:, 1));
raw_eta = squeeze(raw_eta_all(:, :, x_idx));
raw_phi = squeeze(raw_phi_all(:, :, x_idx));
raw_u = squeeze(raw_u_all(:, :, x_idx));
raw_w = squeeze(raw_w_all(:, :, x_idx));

ow3d = struct();
ow3d.time = ref.t(:);
ow3d.eta = four_phase_temporal_separation(raw_eta);
ow3d.phi = four_phase_temporal_separation(raw_phi);
ow3d.u = four_phase_temporal_separation(raw_u);
ow3d.w = four_phase_temporal_separation(raw_w);
ow3d.x_vec = ref.x(:, 1);
ow3d.sigma = sigma_vec;
ow3d.probe_x_index = x_idx;
ow3d.probe_x = ref.x(x_idx, 1);
ow3d.probe_sigma_index = sigma_idx;
ow3d.probe_sigma = sigma_vec(sigma_idx);
ow3d.probe_info = probe_info;
ow3d.case_cfg = case_cfg;
end

function [x_idx, info] = resolve_probe_x_index_local(case_cfg, raw_eta_all, x_vec)
Nt = size(raw_eta_all, 1);
Nx = size(raw_eta_all, 3);
x_vec = x_vec(:);

switch lower(case_cfg.probe_x_mode)
    case 'fixed_index'
        if ~isfield(case_cfg, 'probe_x_index') || isempty(case_cfg.probe_x_index)
            error('case_cfg.probe_x_index is required when probe_x_mode=''fixed_index''.');
        end
        x_idx = min(max(1, round(case_cfg.probe_x_index)), Nx);
        info = struct('mode', 'fixed_index');
    case 'auto_centered'
        center_idx_time = round((Nt + 1) / 2);
        centeredness_weight = 0.35;
        best_score = -Inf;
        x_idx = 1;
        best_peak_idx = 1;
        best_peak_amp = 0;

        for ix = 1:Nx
            eta_sep = four_phase_temporal_separation(squeeze(raw_eta_all(:, :, ix)));
            env = abs(hilbert(eta_sep.first(:)));
            [peak_amp, peak_idx] = max(env);
            centeredness = 1 - abs(peak_idx - center_idx_time) / max(center_idx_time - 1, 1);
            score = peak_amp * (1 + centeredness_weight * centeredness);
            if score > best_score
                best_score = score;
                x_idx = ix;
                best_peak_idx = peak_idx;
                best_peak_amp = peak_amp;
            end
        end

        info = struct( ...
            'mode', 'auto_centered', ...
            'time_center_index', center_idx_time, ...
            'selected_peak_time_index', best_peak_idx, ...
            'selected_peak_amplitude', best_peak_amp, ...
            'selected_center_distance', best_peak_idx - center_idx_time, ...
            'centeredness_weight', centeredness_weight);
    otherwise
        error('Unsupported probe_x_mode: %s', case_cfg.probe_x_mode);
end

info.selected_x_index = x_idx;
info.selected_x = x_vec(x_idx);
end
