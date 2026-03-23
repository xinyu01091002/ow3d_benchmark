%% Compare original dalzell_2d against cleaned local variant

clc; clear; close all;

project_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(project_dir, 'test functions for VWA Opensource'));
addpath('C:\Users\spet5947\OneDrive - Nexus365\Research\code transfer\function');

S = load(fullfile(project_dir, 'processed_boundkinematics', 'OW3D_boundkinematics_raw_tidx_0066.mat'));
phase_names = sort(fieldnames(S.raw_snapshot));
x = S.raw_meta.x(:);
h = mean(S.raw_snapshot.(phase_names{1}).h(:), 'omitnan');

eta_phases = zeros(numel(phase_names), numel(x));
coef = [
    0.25  0    -0.25  0     0    -0.25  0     0.25;
    0.25 -0.25  0.25 -0.25  0     0      0     0;
    0.25  0    -0.25  0     0     0.25  0    -0.25;
    0.25  0.25  0.25  0.25  0     0      0     0];

for idx = 1:numel(phase_names)
    eta_slice = squeeze(S.raw_snapshot.(phase_names{idx}).eta(1, :));
    eta_phases(idx, :) = eta_slice(:).';
end

eta11 = reconstruct_harmonics_1d_local(eta_phases, coef);
eta11 = eta11(1, :);

[eta22_ref, eta20_ref, phi22_ref, phi20_ref] = dalzell_2d(eta11, x, h);
[eta22_new, eta20_new, phi22_new, phi20_new] = dalzell_2d_clean(eta11, x, h);

fprintf('Original dalzell_2d:\n');
fprintf('  max eta20 = %.6e, max phi20 = %.6e\n', max(abs(eta20_ref)), max(abs(phi20_ref)));
fprintf('  max eta22 = %.6e, max phi22 = %.6e\n', max(abs(eta22_ref)), max(abs(phi22_ref)));
fprintf('Clean dalzell_2d:\n');
fprintf('  max eta20 = %.6e, max phi20 = %.6e\n', max(abs(eta20_new)), max(abs(phi20_new)));
fprintf('  max eta22 = %.6e, max phi22 = %.6e\n', max(abs(eta22_new)), max(abs(phi22_new)));

report_metric('eta20', eta20_ref, eta20_new);
report_metric('phi20', phi20_ref, phi20_new);
report_metric('eta22', eta22_ref, eta22_new);
report_metric('phi22', phi22_ref, phi22_new);

function h = reconstruct_harmonics_1d_local(phase_data, coef)
    analytic_part = hilbert(phase_data.').';
    all_fields = zeros(8, size(phase_data, 2));
    all_fields(1:2:end, :) = real(phase_data);
    all_fields(2:2:end, :) = -imag(analytic_part);
    h = zeros(4, size(phase_data, 2));
    for n = 1:4
        h(n, :) = coef(n, :) * all_fields;
    end
end

function report_metric(name, a, b)
    a = a(:);
    b = b(:);
    rmse = sqrt(mean((a - b).^2));
    if all(abs(a) < eps) || all(abs(b) < eps)
        corr_val = NaN;
    else
        c = corrcoef(a, b);
        corr_val = c(1, 2);
    end
    fprintf('  %s: corr = %.6f, rmse = %.6e\n', name, corr_val, rmse);
end
