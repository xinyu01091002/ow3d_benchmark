function out = vwa_compute(eta11, x, d, g, opts)
%VWA_COMPUTE  Variable Wavenumber Approximation (VWA) bound harmonics (2nd-5th)
%
% Implements "Mode II" structure using x-direction analytic signal and a
% single linear filter kappa_n:
%   eta^(nn)   = real( (eta_plus)^(n-1) .* kappa_eta_plus )
%   phi_s^(nn) = imag( (eta_plus)^(n-1) .* kappa_phi_plus )   (default)
%
% Supports:
%   - 1D: eta11 is [Nx x 1] (or [1 x Nx])
%   - 2D: eta11 is [Nx x Ny]
% using the SAME function. All FFT/analytic/filtering is performed along x
% (dimension 1). y-spacing dy may differ from dx; it is not used under the
% directional VWA assumption k -> kx.
%
% Inputs
%   eta11 : real-valued linear surface elevation, size [Nx x Ny]
%   x     : x-coordinate vector, length Nx (uniform spacing assumed)
%   d     : water depth
%   g     : gravity
%   opts  : struct (optional)
%       .nList            : list of orders, default [2 3 4 5]
%       .analytic_side    : 'neg' or 'pos', default 'neg'
%                           'neg' mimics your current implementation:
%                           zero DC and nonnegative kx, keep negative kx, *2
%       .compute_eta      : true/false, default true
%       .compute_phi_s    : true/false, default true
%       .phi_take         : 'imag' or 'real', default 'imag'
%       .small_kd_min     : minimum |kd| used in kernel evaluation, default 1e-12
%
% Output
%   out.eta{n}    : eta^(nn) field [Nx x Ny]
%   out.phi_s{n}  : phi_s^(nn) field [Nx x Ny]
%   out.meta      : metadata (kx, dx, conventions)
%
% Notes
%   - This is the "directional VWA": kernels are evaluated with k = kx only.
%   - No padding/dealiasing by default (as requested). Use at your own risk
%     if bandwidth is wide or grid kmax is insufficient.

    arguments
        eta11 double
        x double
        d (1,1) double {mustBePositive}
        g (1,1) double {mustBePositive}
        opts.nList double = [2 3 4 5]
        opts.analytic_side char = 'neg'
        opts.compute_eta (1,1) logical = true
        opts.compute_phi_s (1,1) logical = true
        opts.phi_take char = 'imag'
        opts.small_kd_min (1,1) double = 1e-12
    end

    % Ensure eta11 is [Nx x Ny]
    if isvector(eta11)
        eta11 = eta11(:);
    end

    Nx = size(eta11, 1);
    Ny = size(eta11, 2);

    x = x(:);
    if numel(x) ~= Nx
        error('vwa_compute: length(x) must equal size(eta11,1).');
    end

    % Check uniform spacing
    dx = mean(diff(x));
    if any(abs(diff(x) - dx) > 1e-10*max(1,abs(dx)))
        error('vwa_compute: x must be (approximately) uniformly spaced for FFT-based method.');
    end

    % kx grid consistent with fft(...,[],1) ordering (no fftshift needed)
    kx = vwa_kxgrid(Nx, dx); % [Nx x 1], signed

    % Precompute analytic signal of eta11 along x
    eta_plus = vwa_analytic_x(eta11, opts.analytic_side);

    % Prepare outputs
    out = struct();
    out.eta   = cell(1, max(opts.nList));
    out.phi_s = cell(1, max(opts.nList));

    % FFT of eta11 along x
    Eta_hat = fft(eta11, [], 1); % [Nx x Ny]

    for n = opts.nList(:).'
        if ~ismember(n, [2 3 4 5])
            error('vwa_compute: nList can only contain 2,3,4,5.');
        end

        % ===== eta^(nn) =====
        if opts.compute_eta
            % Build G_nn(kx) = B_n(kx*d) * kx^(n-1)  (directional: k -> kx)
            Bn = vwa_G_coeff(n, kx, d, opts.small_kd_min); % uses signed kx via kd = kx*d
            Gnn = (kx.^(n-1)) .* Bn;

            % Remove kx=0 singularities cleanly
            Gnn(~isfinite(Gnn)) = 0;

            % Apply linear filter along x: kappa_eta = ifft( Eta_hat .* Gnn )
            kappa_eta = ifft(Eta_hat .* reshape(Gnn, [], 1), [], 1);

            % Analytic of kappa_eta along x
            kappa_eta_plus = vwa_analytic_x(kappa_eta, opts.analytic_side);

            % Mode-II product
            eta_nn = real( eta_plus.^(n-1) .* kappa_eta_plus );
            out.eta{n} = eta_nn;
        end

        % ===== phi_s^(nn) =====
        if opts.compute_phi_s
            % Use mu_nn(k) evaluated with k = |kx| (mu is defined for k>=0)
            mu = vwa_mu(n, abs(kx), d, g, opts.small_kd_min);
            mu(~isfinite(mu)) = 0;

            % Linear filter: kappa_phi = ifft( Eta_hat .* mu )
            kappa_phi = ifft(Eta_hat .* reshape(mu, [], 1), [], 1);

            % Analytic of kappa_phi along x
            kappa_phi_plus = vwa_analytic_x(kappa_phi, opts.analytic_side);

            % Mode-II product: take imag (default) to match your Im{ } convention
            switch lower(opts.phi_take)
                case 'imag'
                    phi_nn = imag( eta_plus.^(n-1) .* kappa_phi_plus );
                case 'real'
                    phi_nn = real( eta_plus.^(n-1) .* kappa_phi_plus );
                otherwise
                    error('vwa_compute: opts.phi_take must be ''imag'' or ''real''.');
            end
            out.phi_s{n} = phi_nn;
        end
    end

    % Metadata
    out.meta = struct();
    out.meta.Nx = Nx;
    out.meta.Ny = Ny;
    out.meta.dx = dx;
    out.meta.kx = kx;
    out.meta.depth = d;
    out.meta.g = g;
    out.meta.analytic_side = opts.analytic_side;
    out.meta.mode = 'directional_k_to_kx';
    out.meta.note = 'All FFT/analytic/filtering performed along x (dim 1).';
end
