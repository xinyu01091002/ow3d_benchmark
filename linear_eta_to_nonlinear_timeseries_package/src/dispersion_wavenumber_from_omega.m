function k = dispersion_wavenumber_from_omega(omega_abs, h, g)
%DISPERSION_WAVENUMBER_FROM_OMEGA Solve w^2 = g k tanh(k h) for k >= 0.

omega_abs = abs(omega_abs);
k = zeros(size(omega_abs));

for idx = 1:numel(omega_abs)
    w = omega_abs(idx);
    if w <= 1e-12
        k(idx) = 0;
        continue;
    end

    mu = (w^2) * h / g;
    if mu < 1
        x = sqrt(mu);
    else
        x = mu;
    end

    for iter = 1:50
        th = tanh(x);
        f = x * th - mu;
        df = th + x * (1 - th^2);
        x_new = x - f / max(df, 1e-12);
        if abs(x_new - x) < 1e-14
            x = x_new;
            break;
        end
        x = x_new;
    end
    k(idx) = x / h;
end
end
