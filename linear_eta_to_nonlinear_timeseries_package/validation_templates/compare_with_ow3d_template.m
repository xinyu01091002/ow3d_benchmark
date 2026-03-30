function compare_with_ow3d_template(result_bundle_path, ow3d_data_path)
%COMPARE_WITH_OW3D_TEMPLATE Template for later manual OW3D validation.
%
% This script is intentionally lightweight. It does not assume any fixed
% OW3D data format beyond "you can load a time vector and matching fields".
% The intended workflow is:
% 1. Load this package result bundle.
% 2. Load an OW3D-derived reference time series prepared by the user.
% 3. Replace the placeholder blocks below with case-specific extraction.

if nargin < 1 || isempty(result_bundle_path)
    error('Provide the path to result_bundle.mat from this package.');
end

S = load(result_bundle_path);
if ~isfield(S, 'result')
    error('Expected variable ''result'' in %s.', result_bundle_path);
end
result = S.result;

fprintf('Loaded package result from: %s\n', result_bundle_path);
if nargin >= 2 && ~isempty(ow3d_data_path)
    fprintf('User-specified OW3D reference path: %s\n', ow3d_data_path);
else
    fprintf('No OW3D reference path supplied yet. This is a template only.\n');
end

disp('Replace the placeholder section in compare_with_ow3d_template.m with your own OW3D readers.');

% -------------------- Placeholder example --------------------
% Example expected reference variables after user loading:
% ref.time
% ref.eta_total
% ref.u_surface
% ref.w_surface
%
% Then you can interpolate and compare, for example:
%
% ref_eta_interp = interp1(ref.time, ref.eta_total, result.time, 'linear', 'extrap');
% corr_eta = corr(result.total.eta(:), ref_eta_interp(:));
% rmse_eta = sqrt(mean((result.total.eta(:) - ref_eta_interp(:)).^2));
% fprintf('eta correlation = %.6f, RMSE = %.6e\n', corr_eta, rmse_eta);
%
% figure('Color', 'w');
% plot(result.time, result.total.eta, 'k-', 'LineWidth', 1.2); hold on;
% plot(result.time, ref_eta_interp, 'r--', 'LineWidth', 1.2);
% legend('package', 'OW3D');
% xlabel('t (s)');
% ylabel('\eta (m)');
% title('Template: package vs OW3D');
% grid on;
end
