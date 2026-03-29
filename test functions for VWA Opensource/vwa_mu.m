function mu = vwa_mu(n, k, d, g, kd_min)
%VWA_MU  mu_{nn}^{VWA}(k) for n=2..5, with k>=0.
%
% omega = sqrt(g k tanh(kd))
% sigma = tanh(kd), alpha1 = cosh(2kd) (implicitly via cosh terms)

    if nargin < 5
        kd_min = 1e-12;
    end

    k = max(k, 0);
    kd = k .* d;

    % Avoid kd=0 exactly
    kd(kd < kd_min) = kd_min;

    omega = sqrt(g .* k .* tanh(kd));
    omega = max(omega, 1e-12);

    csch = @(x) 1./sinh(x);

    switch n
        case 2
            mu = -(omega/8) .* (4 + 3*cosh(2*kd).*csch(kd).^4);

        case 3
            mu = -(k./(64*omega)) .* ( 8*g*k + omega.^2 .* coth(kd) .* ...
                 (16 + 56*csch(kd).^2 + 32*csch(kd).^4 + 9*csch(kd).^6) );

        case 4
            mu = -(k.^2 .* omega) ./ (6144*(2 + 3*cosh(2*kd))) .* csch(kd).^10 .* ...
                 (408 + 638*cosh(2*kd) + 230*cosh(4*kd) + 171*cosh(6*kd) + ...
                  124*cosh(8*kd) + 43*cosh(10*kd) + 6*cosh(12*kd));

        case 5
            mu = -(g*k.^4) ./ (1572864*(2+3*cosh(2*kd)).*(1+4*cosh(2*kd)).*omega) .* csch(kd).^12 .* ...
                 (184650 + 222141*cosh(2*kd) + 182320*cosh(4*kd) + 43815*cosh(6*kd) + ...
                  72200*cosh(8*kd) + 44773*cosh(10*kd) + 20880*cosh(12*kd) + ...
                  6071*cosh(14*kd) + 750*cosh(16*kd));

        otherwise
            error('vwa_mu: n must be 2,3,4,5.');
    end

    mu(~isfinite(mu)) = 0;
    mu(kd < 0.3) = 0;
end
