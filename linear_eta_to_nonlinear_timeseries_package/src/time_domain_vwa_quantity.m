function out = time_domain_vwa_quantity(eta_linear, t, h, g, quantity_name, orders, opts)
%TIME_DOMAIN_VWA_QUANTITY Time-domain VWA reconstruction from eta(t).

if nargin < 6 || isempty(orders)
    orders = [2, 3];
end
if nargin < 7 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'analytic_side') || isempty(opts.analytic_side)
    opts.analytic_side = 'pos';
end
if ~isfield(opts, 'small_kd_min') || isempty(opts.small_kd_min)
    opts.small_kd_min = 0.3;
end

eta_linear = eta_linear(:);
t = t(:);
N = numel(t);
dt = t(2) - t(1);
omega = 2 * pi * [0:ceil(N / 2) - 1, -floor(N / 2):-1]' / (N * dt);
k_abs = dispersion_wavenumber_from_omega(abs(omega), h, g);

eta_plus = hilbert(eta_linear);
fft_eta = fft(eta_linear);

out = struct();
for order = orders(:).'
    switch lower(quantity_name)
        case 'eta'
            Bn = vwa_G_coeff(order, sign(omega) .* k_abs, h, opts.small_kd_min);
            coeff = (sign(omega) .* k_abs).^(order - 1) .* Bn;
            field = ifft(fft(eta_plus) .* coeff);
            out.(sprintf('order%d', order)) = real((eta_plus .^ (order - 1)) .* field);
        case 'phi'
            if order == 1
                coeff = zeros(size(k_abs));
                zero_mask = (k_abs <= 1e-12);
                kd_safe = max(k_abs .* h, 1e-12);
                omega_abs = sqrt(g .* k_abs .* tanh(k_abs .* h));
                coeff(~zero_mask) = -(omega_abs(~zero_mask) ./ k_abs(~zero_mask)) .* coth(kd_safe(~zero_mask));
                coeff(~isfinite(coeff)) = 0;
                coeff(k_abs .* h < opts.small_kd_min) = 0;
                kappa = ifft(fft_eta .* coeff);
                kappa_plus = hilbert(kappa);
                out.(sprintf('order%d', order)) = imag(kappa_plus);
            else
                mu = vwa_mu(order, k_abs, h, g, opts.small_kd_min);
                field = ifft(fft(eta_plus) .* mu);
                out.(sprintf('order%d', order)) = imag((eta_plus .^ (order - 1)) .* field);
            end
        case {'u', 'w'}
            [coeff, phase_type] = vwa_surface_quantity_coeff(quantity_name, order, k_abs, h, g, opts.small_kd_min);
            kappa = ifft(fft_eta .* coeff);
            kappa_plus = hilbert(kappa);

            switch order
                case 1
                    field = kappa_plus;
                case 2
                    field = eta_plus .* kappa_plus;
                case 3
                    field = (eta_plus .^ 2) .* kappa_plus;
                otherwise
                    error('Unsupported order %d for quantity %s.', order, quantity_name);
            end

            out.(sprintf('order%d', order)) = apply_phase_operator_local(field, phase_type);
        otherwise
            error('Unsupported quantity name: %s', quantity_name);
    end
end
end
