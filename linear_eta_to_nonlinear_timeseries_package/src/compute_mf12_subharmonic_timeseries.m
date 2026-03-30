function sub = compute_mf12_subharmonic_timeseries(eta_linear, t, h, g, cfg)
%COMPUTE_MF12_SUBHARMONIC_TIMESERIES Compute MF12 2nd-order difference-frequency terms from eta(t).

if nargin < 5 || isempty(cfg)
    cfg = default_config();
end

spec = build_linear_components_from_time_series(eta_linear, t, h, g, cfg.mf12);
coeffs2 = mf12_direct_coefficients(2, g, h, spec.a.', spec.b.', spec.kx.', spec.ky.', 0, 0, 0);

[eta20, phi20, u20, ~, w20] = mf12_second_subharmonic_kinematics(coeffs2, 0, 0, 0, t(:).');

sub = struct();
sub.eta = eta20(:);
sub.u_surface = u20(:);
sub.w_surface = w20(:);
sub.phi_surface = phi20(:);
sub.meta = spec;
end
