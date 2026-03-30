function out = four_phase_temporal_separation(raw_phases)
%FOUR_PHASE_TEMPORAL_SEPARATION Reconstruct 1st/2nd+/3rd/2- time series.

raw_phases = double(raw_phases);
if size(raw_phases, 2) ~= 4
    error('raw_phases must be [Nt x 4].');
end

four_phase_coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

all_time = [raw_phases, imag(hilbert(raw_phases))];
harmonics = zeros(size(raw_phases, 1), 4);
for idx = 1:4
    harmonics(:, idx) = sum(all_time .* four_phase_coef(idx, :), 2);
end

out = struct();
out.first = harmonics(:, 1);
out.second_super = harmonics(:, 2);
out.third = harmonics(:, 3);
out.second_sub = harmonics(:, 4);
out.total = sum(harmonics, 2);
out.coefficients = four_phase_coef;
end
