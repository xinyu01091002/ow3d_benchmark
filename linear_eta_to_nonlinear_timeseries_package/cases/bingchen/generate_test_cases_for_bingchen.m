function manifest = generate_test_cases_for_bingchen(output_dir)
%GENERATE_TEST_CASES_FOR_BINGCHEN Generate Bingchen-requested wavemaker-style inputs.
%
% The test-case shape follows the shared wavemaker visualization script:
% eta(x,t) = clip(t/4,0,1) * A * cos(k x - omega t) * exp(-((t-t0)/tau)^2)
% and here we store the single-point x = 0 time series eta(t).

if nargin < 1 || isempty(output_dir)
    this_dir = fileparts(mfilename('fullpath'));
    root_dir = fileparts(fileparts(this_dir));
    output_dir = fullfile(root_dir, 'examples', 'bingchen_test_cases');
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

g = 9.81;
h = 10.0;
t0 = 30.0;
tau = 15.0;
dt = 0.1;
tmax = 90.0;
t = (0:dt:tmax).';

case_defs = struct( ...
    'name', {'case_T6_a035_k0047', 'case_T10_a065_k0068', 'case_T14_a100_k0130'}, ...
    'Tp',   {6.0, 10.0, 14.0}, ...
    'A',    {0.35, 0.65, 1.0}, ...
    'kp',   {0.047, 0.068, 0.13});

manifest = repmat(struct(), numel(case_defs), 1);
for idx = 1:numel(case_defs)
    omega = 2 * pi / case_defs(idx).Tp;
    ramp = min(max(t / 4.0, 0.0), 1.0);
    targ = (t - t0) / tau;
    eta_linear = ramp .* case_defs(idx).A .* cos(omega * t) .* exp(-(targ .^ 2));

    output_path = fullfile(output_dir, [case_defs(idx).name '.mat']);
    case_meta = struct( ...
        'source', 'bingchen_wavemaker_visual_shared', ...
        'Tp', case_defs(idx).Tp, ...
        'A', case_defs(idx).A, ...
        'kp', case_defs(idx).kp, ...
        'h', h, ...
        't0', t0, ...
        'tau', tau, ...
        'dt', dt, ...
        'tmax', tmax, ...
        'formula', 'clip(t/4,0,1)*A*cos(omega*t)*exp(-((t-t0)/tau)^2) at x=0');

    kp = case_defs(idx).kp; %#ok<NASGU>
    save(output_path, 't', 'eta_linear', 'h', 'kp', 'g', 'case_meta');

    manifest(idx).name = case_defs(idx).name;
    manifest(idx).input_path = output_path;
    manifest(idx).Tp = case_defs(idx).Tp;
    manifest(idx).A = case_defs(idx).A;
    manifest(idx).kp = case_defs(idx).kp;
    manifest(idx).h = h;
end

manifest_table = struct2table(manifest);
writetable(manifest_table, fullfile(output_dir, 'bingchen_case_manifest.csv'));
fprintf('Bingchen test cases written to:\n%s\n', output_dir);
end
