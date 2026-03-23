function out = vwa_compute_surface_quantity(eta11, x, d, g, quantity_name, opts)
%VWA_COMPUTE_SURFACE_QUANTITY  Directional-VWA surface u/w from linear eta11.
%
% Implements the same x-direction analytic/filter/product structure used by
% the current VWA helper family, but specialized to surface velocity
% quantities supported by the OW3D postprocessing workflow.
%
% Inputs
%   eta11        : linear surface elevation [Nx x 1] or [Nx x Ny]
%   x            : x-coordinate vector, length Nx
%   d            : water depth
%   g            : gravity
%   quantity_name: 'u' or 'w'
%   opts         : options structure
%
% Options
%   .analytic_side : 'neg' or 'pos', default 'neg'
%   .small_kd_min  : minimum kd cutoff, default 0.3

    if nargin < 6 || isempty(opts)
        opts = struct();
    end

    if ~isfield(opts, 'analytic_side') || isempty(opts.analytic_side)
        opts.analytic_side = 'neg';
    end
    if ~isfield(opts, 'small_kd_min') || isempty(opts.small_kd_min)
        opts.small_kd_min = 0.3;
    end

    if ~(isnumeric(d) && isscalar(d) && isfinite(d) && d > 0)
        error('vwa_compute_surface_quantity: d must be a positive scalar.');
    end
    if ~(isnumeric(g) && isscalar(g) && isfinite(g) && g > 0)
        error('vwa_compute_surface_quantity: g must be a positive scalar.');
    end

    if isvector(eta11)
        eta11 = eta11(:);
    end

    Nx = size(eta11, 1);
    Ny = size(eta11, 2);

    x = x(:);
    if numel(x) ~= Nx
        error('vwa_compute_surface_quantity: length(x) must equal size(eta11,1).');
    end

    dx = mean(diff(x));
    if any(abs(diff(x) - dx) > 1e-10 * max(1, abs(dx)))
        error('vwa_compute_surface_quantity: x must be approximately uniformly spaced.');
    end

    kx = vwa_kxgrid(Nx, dx);
    eta_plus = vwa_analytic_x(eta11, opts.analytic_side);
    eta_hat = fft(eta11, [], 1);

    out = struct();
    for order = 1:3
        [coeff, phase_type] = vwa_surface_quantity_coeff(quantity_name, order, abs(kx), d, g, opts.small_kd_min);
        kappa = ifft(eta_hat .* reshape(coeff, [], 1), [], 1);
        kappa_plus = vwa_analytic_x(kappa, opts.analytic_side);

        switch order
            case 1
                field = kappa_plus;
            case 2
                field = eta_plus .* kappa_plus;
            case 3
                field = (eta_plus .^ 2) .* kappa_plus;
        end

        out.(sprintf('order%d', order)) = vwa_apply_phase_operator(field, phase_type);
    end

    out.meta = struct();
    out.meta.Nx = Nx;
    out.meta.Ny = Ny;
    out.meta.dx = dx;
    out.meta.kx = kx;
    out.meta.depth = d;
    out.meta.g = g;
    out.meta.quantity_name = lower(quantity_name);
    out.meta.analytic_side = opts.analytic_side;
    out.meta.small_kd_min = opts.small_kd_min;
    out.meta.mode = 'directional_k_to_kx_surface_quantity';
end
