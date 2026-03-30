root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'src'));
package_setup();
fprintf('linear_eta_to_nonlinear_timeseries_package is now on the MATLAB path.\n');
