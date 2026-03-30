function [coeff, phase_type] = vwa_surface_quantity_coeff(quantity_name, order, k_abs, d, g, kd_min)
%VWA_SURFACE_QUANTITY_COEFF  Directional-VWA transfer coefficients for surface quantities.
%
% Currently supported quantities:
%   - 'u' : surface horizontal velocity
%   - 'w' : surface vertical velocity
%
% The formulas are extracted from the current OW3D postprocessing workflow
% so the reusable VWA helpers match the existing in-script implementation.

    if nargin < 6
        kd_min = 0.3;
    end

    kd = k_abs .* d;
    kd_safe = max(kd, 1e-12);
    sigma = tanh(kd_safe);
    omega = sqrt(g .* k_abs .* sigma);

    coeff = zeros(size(k_abs));
    zero_mask = (k_abs <= 1e-12);
    phase_type = 'real';

    switch lower(quantity_name)
        case 'u'
            phase_type = 'real';
            switch order
                case 1
                    coeff(~zero_mask) = omega(~zero_mask) .* coth(kd_safe(~zero_mask));
                case 2
                    coeff(~zero_mask) = k_abs(~zero_mask) .* omega(~zero_mask) .* ...
                        (0.5 + 3 * cosh(2 * kd_safe(~zero_mask)) ./ (4 * sinh(kd_safe(~zero_mask)).^4));
                case 3
                    term1 = -(3 / 64) * k_abs.^2 .* omega .* ...
                        ((-11 + 2 * cosh(2 * kd_safe)) ./ sinh(kd_safe).^7) .* cosh(3 * kd_safe);
                    term2 = (3 / 2) * k_abs.^2 .* omega .* cosh(kd_safe) ./ sinh(kd_safe).^3;
                    term3 = (1 / 8) * k_abs.^2 .* omega .* ...
                        (2 + cosh(2 * kd_safe)) .* cosh(kd_safe) ./ sinh(kd_safe).^3;
                    term4 = (1 / 8) * omega .* k_abs.^2 .* coth(kd_safe);
                    coeff(~zero_mask) = term1(~zero_mask) + term2(~zero_mask) + term3(~zero_mask) + term4(~zero_mask);
                otherwise
                    error('vwa_surface_quantity_coeff: unsupported order %d for quantity %s.', order, quantity_name);
            end

        case 'w'
            phase_type = 'imag';
            switch order
                case 1
                    coeff(~zero_mask) = -omega(~zero_mask);
                case 2
                    coeff(~zero_mask) = -omega(~zero_mask) .* k_abs(~zero_mask) .* ...
                        (0.5 .* coth(kd_safe(~zero_mask)) + 1.5 .* cosh(kd_safe(~zero_mask)) ./ sinh(kd_safe(~zero_mask)).^3);
                case 3
                    term1 = (3 / 64) * omega .* k_abs.^2 .* ...
                        ((-11 + 2 * cosh(2 * kd_safe)) ./ sinh(kd_safe).^7) .* sinh(3 * kd_safe);
                    term2 = -(3 / 4) * omega .* k_abs.^2 .* cosh(2 * kd_safe) ./ sinh(kd_safe).^4;
                    term3 = -(1 / 8) * omega .* k_abs.^2 .* ...
                        (2 + cosh(2 * kd_safe)) .* coth(kd_safe).^2 ./ sinh(kd_safe).^2;
                    term4 = -(1 / 8) * omega .* k_abs.^2;
                    coeff(~zero_mask) = term1(~zero_mask) + term2(~zero_mask) + term3(~zero_mask) + term4(~zero_mask);
                otherwise
                    error('vwa_surface_quantity_coeff: unsupported order %d for quantity %s.', order, quantity_name);
            end

        otherwise
            error('vwa_surface_quantity_coeff: unsupported quantity ''%s''.', quantity_name);
    end

    coeff(~isfinite(coeff)) = 0;
    coeff(zero_mask) = 0;
    coeff(kd < kd_min) = 0;
end
