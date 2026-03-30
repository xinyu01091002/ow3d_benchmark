function ow3d = extract_ow3d_unidirectional_timeseries_case(case_cfg)
%EXTRACT_OW3D_UNIDIRECTIONAL_TIMESERIES_CASE Extract a four-phase OW3D time series case.

required = {'data_root', 'Alpha', 'Akp', 'kd', 't_steps', 'dt', 'target_idx'};
for k = 1:numel(required)
    if ~isfield(case_cfg, required{k})
        error('case_cfg must contain field ''%s''.', required{k});
    end
end

phi_shifts = 0:90:270;
Nt = numel(case_cfg.t_steps);
raw_eta = zeros(Nt, 4);
raw_phi = zeros(Nt, 4);
x_vec = [];
y_vec = [];

for p_idx = 1:numel(phi_shifts)
    phi_shift = phi_shifts(p_idx);
    folder_name = sprintf('T_init-40_Tp_Alpha_%.1f_Akp_%.3d_kd%.1f_phi_%d', ...
        case_cfg.Alpha, round(100 * case_cfg.Akp), case_cfg.kd, phi_shift);
    folder_path = fullfile(case_cfg.data_root, folder_name);
    if ~isfolder(folder_path)
        error('Missing OW3D folder: %s', folder_path);
    end

    for t_idx = 1:Nt
        step = case_cfg.t_steps(t_idx);
        file_path = fullfile(folder_path, sprintf('EP_%05d.bin', step));
        if ~isfile(file_path)
            error('Missing OW3D snapshot: %s', file_path);
        end

        [X, Y, eta_tmp, phi_tmp] = read_ow3d_surface_bin(file_path);
        if isempty(x_vec)
            x_vec = X(:, 1);
            y_vec = Y(1, :);
        end
        raw_eta(t_idx, p_idx) = eta_tmp(case_cfg.target_idx);
        raw_phi(t_idx, p_idx) = phi_tmp(case_cfg.target_idx);
    end
end

eta_sep = four_phase_temporal_separation(raw_eta);
phi_sep = four_phase_temporal_separation(raw_phi);

ow3d = struct();
ow3d.time = case_cfg.t_steps(:) * case_cfg.dt;
ow3d.raw_eta_phases = raw_eta;
ow3d.raw_phi_phases = raw_phi;
ow3d.eta = eta_sep;
ow3d.phi = phi_sep;
ow3d.x_vec = x_vec;
ow3d.y_vec = y_vec;
ow3d.target_idx = case_cfg.target_idx;
ow3d.target_x = x_vec(case_cfg.target_idx);
ow3d.case_cfg = case_cfg;
end
