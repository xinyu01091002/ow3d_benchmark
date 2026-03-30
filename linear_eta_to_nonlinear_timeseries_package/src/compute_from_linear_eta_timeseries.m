function result = compute_from_linear_eta_timeseries(input_data, cfg)
%COMPUTE_FROM_LINEAR_ETA_TIMESERIES Top-level package API.

if nargin < 2 || isempty(cfg)
    cfg = default_config();
end

if ischar(input_data) || isstring(input_data)
    input_data = load_linear_eta_input(input_data, cfg);
end

t = input_data.t(:);
eta_linear = input_data.eta_linear(:);
h = input_data.h;
kp = input_data.kp;
g = input_data.g;

linear = compute_linear_surface_velocity_timeseries(eta_linear, t, h, g, cfg);
super = compute_vwa_superharmonic_timeseries(eta_linear, t, h, g, cfg);
sub = compute_mf12_subharmonic_timeseries(eta_linear, t, h, g, cfg);
total = assemble_total_timeseries(linear, super, sub);

result = struct();
result.time = t;
result.linear = struct( ...
    'eta', linear.eta(:), ...
    'phi_surface', linear.phi_surface(:), ...
    'u_surface', linear.u_surface(:), ...
    'w_surface', linear.w_surface(:));
result.superharmonic = struct( ...
    'eta', super.eta(:), ...
    'phi_surface', super.phi_surface(:), ...
    'u_surface', super.u_surface(:), ...
    'w_surface', super.w_surface(:), ...
    'orders', super.orders);
result.subharmonic = struct( ...
    'eta', sub.eta(:), ...
    'phi_surface', sub.phi_surface(:), ...
    'u_surface', sub.u_surface(:), ...
    'w_surface', sub.w_surface(:));
result.total = total;
result.nonlinear = total;
result.meta = struct();
result.meta.source_file = input_data.source_file;
result.meta.h = h;
result.meta.kp = kp;
result.meta.g = g;
result.meta.dt = input_data.dt;
result.meta.num_samples = numel(t);
result.meta.Tp = 2 * pi / sqrt(g * kp * tanh(kp * h));
result.meta.package_root = fileparts(fileparts(mfilename('fullpath')));
result.meta.timestamp = datestr(now, 0);
result.meta.assumptions = [ ...
    "Input is a unidirectional single-point linear eta(t)."; ...
    "VWA superharmonics are computed in time domain using the existing coefficient family."; ...
    "MF12 subharmonics are computed from a frequency decomposition of eta(t) at a fixed spatial point x=0."; ...
    "Surface velocities are exported as the package's default physical-surface outputs, but still require case-by-case scientific validation."];
result.debug = struct('linear', linear, 'superharmonic', super, 'subharmonic', sub);
end
