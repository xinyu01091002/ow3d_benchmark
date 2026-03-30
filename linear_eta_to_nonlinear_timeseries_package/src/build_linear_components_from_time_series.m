function linear_spec = build_linear_components_from_time_series(eta_linear, t, h, g, opts)
%BUILD_LINEAR_COMPONENTS_FROM_TIME_SERIES Convert eta(t) to MF12-style linear components.

if nargin < 5 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'energy_keep') || isempty(opts.energy_keep)
    opts.energy_keep = 1.0;
end
if ~isfield(opts, 'max_components') || isempty(opts.max_components)
    opts.max_components = inf;
end
if ~isfield(opts, 'drop_zero_frequency') || isempty(opts.drop_zero_frequency)
    opts.drop_zero_frequency = true;
end

eta_linear = eta_linear(:);
t = t(:);
N = numel(eta_linear);
dt = t(2) - t(1);

fft_eta = fft(eta_linear) / N;
omega = 2 * pi * [0:ceil(N / 2) - 1, -floor(N / 2):-1]' / (N * dt);

if mod(N, 2) == 0
    positive_idx = 2:(N / 2);
else
    positive_idx = 2:((N + 1) / 2);
end

amps = abs(fft_eta(positive_idx));
nonzero = amps > 1e-12 * max(amps + eps);
positive_idx = positive_idx(nonzero);

if opts.drop_zero_frequency
    positive_idx = positive_idx(abs(omega(positive_idx)) > 1e-12);
end

if isempty(positive_idx)
    error('No usable positive-frequency components were found in the input eta(t).');
end

energy = abs(fft_eta(positive_idx)).^2;
total_energy = sum(energy);
if total_energy > 0 && opts.energy_keep < 1
    [energy_sorted, order] = sort(energy, 'descend');
    cutoff = find(cumsum(energy_sorted) >= opts.energy_keep * total_energy, 1, 'first');
    keep = positive_idx(order(1:cutoff));
else
    keep = positive_idx;
end

if isfinite(opts.max_components)
    energy_keep = abs(fft_eta(keep)).^2;
    [~, order] = sort(energy_keep, 'descend');
    keep = keep(order(1:min(numel(order), opts.max_components)));
end

keep = sort(keep(:));
coeff_pos = fft_eta(keep);
omega_pos = omega(keep);
k_pos = dispersion_wavenumber_from_omega(abs(omega_pos), h, g);

linear_spec = struct();
linear_spec.indices = keep(:);
linear_spec.fft_coeff = coeff_pos(:);
linear_spec.omega = abs(omega_pos(:));
linear_spec.k = k_pos(:);
linear_spec.a = 2 * real(coeff_pos(:));
linear_spec.b = -2 * imag(coeff_pos(:));
linear_spec.kx = k_pos(:);
linear_spec.ky = zeros(size(k_pos(:)));
linear_spec.energy_keep = opts.energy_keep;
linear_spec.max_components = opts.max_components;
linear_spec.dt = dt;
linear_spec.num_components = numel(keep);
end
