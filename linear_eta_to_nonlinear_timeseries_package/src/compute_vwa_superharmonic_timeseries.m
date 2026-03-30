function super = compute_vwa_superharmonic_timeseries(eta_linear, t, h, g, cfg)
%COMPUTE_VWA_SUPERHARMONIC_TIMESERIES Compute 2nd/3rd-order superharmonics from eta(t).

if nargin < 5 || isempty(cfg)
    cfg = default_config();
end

orders = cfg.vwa.superharmonic_orders;
opts = struct( ...
    'analytic_side', cfg.vwa.analytic_side, ...
    'small_kd_min', cfg.vwa.small_kd_min);

eta_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'eta', orders, opts);
phi_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'phi', orders, opts);
u_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'u', orders, opts);
w_terms = time_domain_vwa_quantity(eta_linear, t, h, g, 'w', orders, opts);

super = struct();
super.eta = zeros(numel(t), 1);
super.phi_surface = zeros(numel(t), 1);
super.u_surface = zeros(numel(t), 1);
super.w_surface = zeros(numel(t), 1);
super.orders = struct();

for order = orders(:).'
    tag = sprintf('order%d', order);
    super.orders.(tag) = struct( ...
        'eta', eta_terms.(tag)(:), ...
        'phi_surface', phi_terms.(tag)(:), ...
        'u_surface', u_terms.(tag)(:), ...
        'w_surface', w_terms.(tag)(:));
    super.eta = super.eta + super.orders.(tag).eta;
    super.phi_surface = super.phi_surface + super.orders.(tag).phi_surface;
    super.u_surface = super.u_surface + super.orders.(tag).u_surface;
    super.w_surface = super.w_surface + super.orders.(tag).w_surface;
end

super.meta = struct('orders', orders(:).', 'method', 'time_domain_vwa');
end
