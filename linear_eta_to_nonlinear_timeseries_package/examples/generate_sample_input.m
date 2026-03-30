function output_path = generate_sample_input(output_path)
%GENERATE_SAMPLE_INPUT Generate a small demo input .mat file.

if nargin < 1 || isempty(output_path)
    output_path = fullfile(fileparts(mfilename('fullpath')), 'sample_linear_eta_input.mat');
end

g = 9.81;
kp = 0.0279;
h = 5 / kp;
wp = sqrt(g * kp * tanh(kp * h));
Tp = 2 * pi / wp;

dt = Tp / 60;
t = (-20 * Tp:dt:20 * Tp).';
envelope = exp(-(t / (4.0 * Tp)).^2);
eta_linear = 0.75 * (0.02 / kp) * envelope .* cos(wp * t + 0.15) ...
    + 0.18 * (0.02 / kp) * envelope .* cos(1.12 * wp * t - 0.60);

save(output_path, 't', 'eta_linear', 'h', 'kp', 'g');
fprintf('Sample input written to: %s\n', output_path);
end
