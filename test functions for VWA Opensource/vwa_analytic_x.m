function xA = vwa_analytic_x(x, side)
%VWA_ANALYTIC_X  "Analytic" complex field along x (dimension 1) via FFT masking.
%
% side='neg' (default in vwa_compute):
%   - Zero DC and nonnegative kx bins
%   - Keep negative-kx bins, multiply by 2
%
% side='pos':
%   - Zero DC and nonpositive kx bins
%   - Keep positive-kx bins, multiply by 2
%
% This is a pragmatic one-sided spectrum construction used in your current scripts.
% It is not the same as MATLAB's hilbert() "standard" analytic signal, because
% DC (and Nyquist when present) are set to zero here.

    if nargin < 2
        side = 'neg';
    end

    Nx = size(x, 1);

    X = fft(x, [], 1);

    mask = zeros(Nx, 1);

    switch lower(side)
        case 'neg'
            % Keep negative frequencies:
            % indices: floor(N/2)+2 : N
            i0 = floor(Nx/2) + 2;
            if i0 <= Nx
                mask(i0:Nx) = 2;
            end

        case 'pos'
            % Keep positive frequencies:
            % indices: 2 : floor(N/2)+1
            i1 = 2;
            i2 = floor(Nx/2) + 1;
            if i1 <= i2
                mask(i1:i2) = 2;
            end

        otherwise
            error('vwa_analytic_x: side must be ''neg'' or ''pos''.');
    end

    xA = ifft(X .* reshape(mask, [], 1), [], 1);
end
