function run_bingchen_test_case_batch()
%RUN_BINGCHEN_TEST_CASE_BATCH Generate and run Bingchen-requested test cases.

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
addpath(this_dir);
addpath(root_dir);
addpath(fullfile(root_dir, 'src'));
package_setup();

input_dir = fullfile(root_dir, 'examples', 'bingchen_test_cases');
output_root = fullfile(root_dir, 'output_example', 'bingchen_test_cases');
manifest = generate_test_cases_for_bingchen(input_dir);

if ~isfolder(output_root)
    mkdir(output_root);
end

summary_rows = strings(numel(manifest), 5);
for idx = 1:numel(manifest)
    case_output_dir = fullfile(output_root, manifest(idx).name);
    result = run_from_user_eta_timeseries(manifest(idx).input_path, case_output_dir); %#ok<NASGU>
    summary_rows(idx, :) = [ ...
        string(manifest(idx).name), ...
        string(manifest(idx).Tp), ...
        string(manifest(idx).A), ...
        string(manifest(idx).kp), ...
        string(manifest(idx).h)];
end

summary_path = fullfile(output_root, 'bingchen_case_summary.csv');
fid = fopen(summary_path, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'name,Tp,A,kp,h\n');
for idx = 1:size(summary_rows, 1)
    fprintf(fid, '%s,%s,%s,%s,%s\n', summary_rows(idx, 1), summary_rows(idx, 2), summary_rows(idx, 3), summary_rows(idx, 4), summary_rows(idx, 5));
end

plot_bingchen_test_case_overview();

fprintf('Bingchen batch completed. Output written to:\n%s\n', output_root);
end
