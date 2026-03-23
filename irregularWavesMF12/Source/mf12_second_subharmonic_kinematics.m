function [eta, phi, u, v, w, p] = mf12_second_subharmonic_kinematics(coeffs, x, y, z, t)
%MF12_SECOND_SUBHARMONIC_KINEMATICS Evaluate MF12 second-order difference terms only.
%
% This helper extracts the second-order subharmonic (difference-frequency)
% part of the MF12 solution, including free-surface elevation, velocity
% potential, velocities, and dynamic pressure.
%
% Inputs
%   coeffs : structure returned by mf12_direct_coefficients or
%            mf12_spectral_coefficients with order >= 2
%   x,y,z,t: evaluation coordinates/scalars/arrays compatible with MF12
%
% Outputs
%   eta : second-order difference-wave free-surface elevation
%   phi : second-order difference-wave velocity potential
%   u,v,w : velocity components from difference-wave potential
%   p : dynamic pressure p^+/rho from difference-wave potential

    if ~isfield(coeffs, 'omega_npm') || ~isfield(coeffs, 'F_npm') || ~isfield(coeffs, 'G_npm')
        error('mf12_second_subharmonic_kinematics: coeffs must contain second-order pair-interaction fields.');
    end

    eta = zeros(size(t));
    phi = zeros(size(t));
    u = zeros(size(t));
    v = zeros(size(t));
    w = zeros(size(t));
    p = zeros(size(t));

    Z = z + coeffs.h;
    cnm = 0;
    for n = 1:coeffs.N
        for m = (n + 1):coeffs.N
            % In MF12 coefficient arrays, odd entry = n+m, even entry = n-m.
            cnm = cnm + 1;
            idxMinus = 2 * cnm;

            omega_nm = coeffs.omega_npm(idxMinus);
            kx_nm = coeffs.kx_npm(idxMinus);
            ky_nm = coeffs.ky_npm(idxMinus);
            kappa_nm = coeffs.kappa_npm(idxMinus);
            F_nm = coeffs.F_npm(idxMinus);
            G_nm = coeffs.G_npm(idxMinus);
            A_nm = coeffs.A_npm(idxMinus);
            B_nm = coeffs.B_npm(idxMinus);
            mu_nm = coeffs.mu_npm(idxMinus);

            theta_nm = omega_nm .* t - kx_nm .* x - ky_nm .* y;

            eta = eta + G_nm .* (A_nm .* cos(theta_nm) + B_nm .* sin(theta_nm));

            factorZ = F_nm .* cosh(kappa_nm .* Z);
            phi = phi + factorZ .* (A_nm .* sin(theta_nm) - B_nm .* cos(theta_nm));

            u = u + kx_nm .* factorZ .* (-A_nm .* cos(theta_nm) - B_nm .* sin(theta_nm));
            v = v + ky_nm .* factorZ .* (-A_nm .* cos(theta_nm) - B_nm .* sin(theta_nm));
            w = w + F_nm .* sinh(kappa_nm .* Z) .* kappa_nm .* (A_nm .* sin(theta_nm) - B_nm .* cos(theta_nm));
            p = p - factorZ .* omega_nm .* (A_nm .* cos(theta_nm) + B_nm .* sin(theta_nm));
        end
    end

    % Free-surface potential consistent with mf12_direct_surface uses mu_npm.
    if all(abs(Z - coeffs.h) < 1e-12, 'all')
        phi = zeros(size(t));
        cnm = 0;
        for n = 1:coeffs.N
            for m = (n + 1):coeffs.N
                cnm = cnm + 1;
                idxMinus = 2 * cnm;
                theta_nm = coeffs.omega_npm(idxMinus) .* t ...
                    - coeffs.kx_npm(idxMinus) .* x - coeffs.ky_npm(idxMinus) .* y;
                phi = phi + coeffs.mu_npm(idxMinus) .* ...
                    (coeffs.A_npm(idxMinus) .* sin(theta_nm) - coeffs.B_npm(idxMinus) .* cos(theta_nm));
            end
        end
    end
end
