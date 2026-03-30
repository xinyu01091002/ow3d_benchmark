function input_data = load_linear_eta_input(input_path, cfg)
%LOAD_LINEAR_ETA_INPUT Load and validate a user .mat input file.

if nargin < 2 || isempty(cfg)
    cfg = default_config();
end

if ~(ischar(input_path) || isstring(input_path))
    error('Input path must be a character vector or string.');
end
if ~isfile(input_path)
    error('Input file not found: %s', input_path);
end

S = load(input_path);
required = {'t', 'eta_linear', 'h', 'kp'};
for k = 1:numel(required)
    if ~isfield(S, required{k})
        error('Input file must contain field ''%s''.', required{k});
    end
end

input_data = struct();
input_data.t = S.t(:);
input_data.eta_linear = S.eta_linear(:);
input_data.h = double(S.h);
input_data.kp = double(S.kp);

if isfield(S, 'g') && ~isempty(S.g)
    input_data.g = double(S.g);
else
    input_data.g = cfg.g;
end

if numel(input_data.t) ~= numel(input_data.eta_linear)
    error('Fields ''t'' and ''eta_linear'' must have the same length.');
end
if ~isscalar(input_data.h) || ~isfinite(input_data.h) || input_data.h <= 0
    error('Field ''h'' must be a positive scalar.');
end
if ~isscalar(input_data.kp) || ~isfinite(input_data.kp) || input_data.kp <= 0
    error('Field ''kp'' must be a positive scalar.');
end
if ~isscalar(input_data.g) || ~isfinite(input_data.g) || input_data.g <= 0
    error('Field ''g'' must be a positive scalar.');
end
if any(~isfinite(input_data.t)) || any(~isfinite(input_data.eta_linear))
    error('Input time series must not contain NaN or Inf values.');
end
if numel(input_data.t) < 16
    error('Input time series is too short. Provide at least 16 samples.');
end

dt = diff(input_data.t);
dt_ref = mean(dt);
if dt_ref <= 0
    error('Time vector must be strictly increasing.');
end
if any(abs(dt - dt_ref) > 1e-8 * max(1, abs(dt_ref)))
    error('Time vector must be approximately uniformly spaced.');
end

input_data.dt = dt_ref;
input_data.fs = 1 / dt_ref;
input_data.source_file = char(input_path);
end
