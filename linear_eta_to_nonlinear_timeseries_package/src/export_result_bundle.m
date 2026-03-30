function export_result_bundle(result, output_dir, cfg)
%EXPORT_RESULT_BUNDLE Save result bundle to MAT/CSV/PNG files.

if nargin < 3 || isempty(cfg)
    cfg = default_config();
end

if ~isfolder(output_dir)
    mkdir(output_dir);
end

save(fullfile(output_dir, 'result_bundle.mat'), 'result');

if cfg.save_csv
    T = table( ...
        result.time(:), ...
        result.linear.eta(:), result.superharmonic.eta(:), result.subharmonic.eta(:), result.nonlinear.eta(:), ...
        result.linear.phi_surface(:), result.superharmonic.phi_surface(:), result.subharmonic.phi_surface(:), result.nonlinear.phi_surface(:), ...
        result.linear.u_surface(:), result.superharmonic.u_surface(:), result.subharmonic.u_surface(:), result.nonlinear.u_surface(:), ...
        result.linear.w_surface(:), result.superharmonic.w_surface(:), result.subharmonic.w_surface(:), result.nonlinear.w_surface(:), ...
        'VariableNames', { ...
            'time', ...
            'eta_linear', 'eta_superharmonic', 'eta_subharmonic', 'eta_nonlinear', ...
            'phi_linear', 'phi_superharmonic', 'phi_subharmonic', 'phi_nonlinear', ...
            'u_linear', 'u_superharmonic', 'u_subharmonic', 'u_nonlinear', ...
            'w_linear', 'w_superharmonic', 'w_subharmonic', 'w_nonlinear'});
    writetable(T, fullfile(output_dir, 'result_timeseries.csv'));
end

figs = plot_result_summary(result, output_dir, cfg); %#ok<NASGU>

if cfg.save_debug
    debug_summary = struct();
    debug_summary.max_abs_eta_closure = max(abs(result.nonlinear.eta - (result.linear.eta + result.superharmonic.eta + result.subharmonic.eta)));
    debug_summary.max_abs_phi_closure = max(abs(result.nonlinear.phi_surface - (result.linear.phi_surface + result.superharmonic.phi_surface + result.subharmonic.phi_surface)));
    debug_summary.max_abs_u_closure = max(abs(result.nonlinear.u_surface - (result.linear.u_surface + result.superharmonic.u_surface + result.subharmonic.u_surface)));
    debug_summary.max_abs_w_closure = max(abs(result.nonlinear.w_surface - (result.linear.w_surface + result.superharmonic.w_surface + result.subharmonic.w_surface)));
    save(fullfile(output_dir, 'debug_summary.mat'), 'debug_summary');
end
end
