function result = run_from_user_eta_timeseries(input_path, output_dir)
%RUN_FROM_USER_ETA_TIMESERIES Main user entry point for a .mat eta(t) input.

if nargin < 1 || isempty(input_path)
    error(['Usage: run_from_user_eta_timeseries(''', ...
        'C:\path\to\your_input.mat'', ''optional_output_dir'')']);
end

root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'src'));
package_setup();

cfg = default_config();
if nargin >= 2 && ~isempty(output_dir)
    cfg.output_dir = output_dir;
else
    [~, base_name] = fileparts(char(input_path));
    cfg.output_dir = fullfile(root_dir, 'output_example', ['run_' base_name]);
end

input_data = load_linear_eta_input(input_path, cfg);
result = compute_from_linear_eta_timeseries(input_data, cfg);
export_result_bundle(result, cfg.output_dir, cfg);

fprintf('Package run completed. Output written to:\n%s\n', cfg.output_dir);
end
