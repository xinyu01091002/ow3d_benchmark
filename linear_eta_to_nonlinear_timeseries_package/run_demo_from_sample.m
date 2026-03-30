function run_demo_from_sample()
%RUN_DEMO_FROM_SAMPLE Generate a demo input and run the full package.

root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'src'));
package_setup();

cfg = default_config();
cfg.output_dir = fullfile(root_dir, 'output_example', 'demo_run');

input_path = fullfile(root_dir, 'examples', 'sample_linear_eta_input.mat');
if ~isfile(input_path)
    generate_sample_input(input_path);
end

input_data = load_linear_eta_input(input_path, cfg);
result = compute_from_linear_eta_timeseries(input_data, cfg); %#ok<NASGU>
export_result_bundle(result, cfg.output_dir, cfg);

fprintf('Demo completed. Output written to:\n%s\n', cfg.output_dir);
end
