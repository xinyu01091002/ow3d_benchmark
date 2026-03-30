function total = assemble_total_timeseries(linear, super, sub)
%ASSEMBLE_TOTAL_TIMESERIES Assemble linear + super + sub results.

total = struct();
total.eta = linear.eta(:) + super.eta(:) + sub.eta(:);
total.phi_surface = linear.phi_surface(:) + super.phi_surface(:) + sub.phi_surface(:);
total.u_surface = linear.u_surface(:) + super.u_surface(:) + sub.u_surface(:);
total.w_surface = linear.w_surface(:) + super.w_surface(:) + sub.w_surface(:);
end
