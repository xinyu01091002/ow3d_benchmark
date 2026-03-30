function input_data = build_input_from_ow3d_linear_timeseries(ow3d, g)
%BUILD_INPUT_FROM_OW3D_LINEAR_TIMESERIES Build package input from extracted OW3D eta^(1).

if nargin < 2 || isempty(g)
    g = 9.81;
end

input_data = struct();
input_data.t = ow3d.time(:);
input_data.eta_linear = ow3d.eta.first(:);
input_data.h = ow3d.case_cfg.kd / ow3d.case_cfg.kp;
input_data.kp = ow3d.case_cfg.kp;
input_data.g = g;
input_data.dt = mean(diff(input_data.t));
input_data.fs = 1 / input_data.dt;
input_data.source_file = sprintf('OW3D:%s', ow3d.case_cfg.data_root);
end
