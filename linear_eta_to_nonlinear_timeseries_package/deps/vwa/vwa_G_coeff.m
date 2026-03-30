function B = vwa_G_coeff(n, kx, d, kd_min)
%VWA_G_COEFF  Dimensionless coefficient B_n(kd) for G_nn^{VWA}.
%
% Implements:
%   G_{2n} = ((3 - sigma^2)/(4 sigma^3)) * k
%   G_{3n} = ((27 - 9 sigma^2 + 9 sigma^4 - 3 sigma^6)/(64 sigma^6)) * k^2
%   G_{4n} = [poly(alpha1)] / [24 (3 alpha1+2) (alpha1-1)^4 sinh(2kd)] * k^3
%   G_{5n} = [poly(alpha1)] / [384 (alpha1-1)^6 (12 alpha1^2 + 11 alpha1 + 2)] * k^4
%
% Here sigma=tanh(kd), alpha1=cosh(2kd).
%
% Important:
%   We evaluate with kd = kx*d (signed). This matches your current practice
%   where sign changes in sinh(2kd) can cancel with kx^(n-1).
%
% Output B is dimensionless; your main code forms Gnn = kx^(n-1) .* B.

    if nargin < 4
        kd_min = 1e-12;
    end

    kd = kx .* d;

    % Clip only extremely small kd to avoid division blow-up (preserve sign)
    small = abs(kd) < kd_min;
    kd(small) = sign(kd(small) + (kd(small)==0)) * kd_min;

    sigma  = tanh(kd);
    alpha1 = cosh(2*kd);

    % Prevent sigma=0
    sigma(abs(sigma) < 1e-12) = 1e-12 .* sign(sigma(abs(sigma) < 1e-12) + 1);

    switch n
        case 2
            B = (3 - sigma.^2) ./ (4 * sigma.^3);

        case 3
            B = (27 - 9*sigma.^2 + 9*sigma.^4 - 3*sigma.^6) ./ (64 * sigma.^6);

        case 4
            num = 24*alpha1.^6 + 116*alpha1.^5 + 214*alpha1.^4 + 188*alpha1.^3 + ...
                  133*alpha1.^2 + 101*alpha1 + 34;
            den = 24 .* (3*alpha1 + 2) .* (alpha1 - 1).^4 .* sinh(2*kd);
            den = sign(den) .* max(abs(den), 1e-12);
            B = num ./ den;

        case 5
            num = 5*(300*alpha1.^8 + 1579*alpha1.^7 + 3176*alpha1.^6 + 2949*alpha1.^5 + ...
                     1188*alpha1.^4 + 675*alpha1.^3 + 1326*alpha1.^2 + 827*alpha1 + 130);
            den = 384 .* (alpha1 - 1).^6 .* (12*alpha1.^2 + 11*alpha1 + 2);
            den = sign(den) .* max(abs(den), 1e-12);
            B = num ./ den;

        otherwise
            error('vwa_G_coeff: n must be 2,3,4,5.');
    end

    B(~isfinite(B)) = 0;
    B(abs(kd) < 0.3) = 0;
end
