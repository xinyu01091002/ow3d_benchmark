function linear = compute_linear_surface_velocity_timeseries(eta_linear, t, h, g, cfg)
%COMPUTE_LINEAR_SURFACE_VELOCITY_TIMESERIES Compute linear eta/u/w.
%
% For this package we use the same time-domain one-sided construction used
% by the VWA bridge for order-1 u/w, rather than rebuilding an MF12
% periodic component set from a finite window eta(t). The MF12 order-1
% reconstruction is mathematically fine for periodic signals, but for the
% finite-window demo/use case here it can push peak u/w activity to the
% record edges.

if nargin < 5 || isempty(cfg)
    cfg = default_config();
end

opts = struct( ...
    'analytic_side', cfg.vwa.analytic_side, ...
    'small_kd_min', cfg.vwa.small_kd_min);
phi_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'phi', 1, opts);
u_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'u', 1, opts);
w_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'w', 1, opts);

linear = struct();
linear.eta = eta_linear(:);
linear.phi_surface = phi_terms.order1(:);
linear.u_surface = u_terms.order1(:);
linear.w_surface = w_terms.order1(:);
linear.meta = struct('method', 'time_domain_order1_relation', 'phi_method', 'time_domain_transfer_coeff');
end
