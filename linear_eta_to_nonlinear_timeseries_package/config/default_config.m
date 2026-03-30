function cfg = default_config()
%DEFAULT_CONFIG Default configuration for the linear-eta package.

cfg = struct();

cfg.g = 9.81;
cfg.output_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'output_example');
cfg.save_csv = true;
cfg.save_figures = true;
cfg.save_debug = true;

cfg.vwa = struct();
cfg.vwa.superharmonic_orders = [2, 3];
cfg.vwa.analytic_side = 'pos';
cfg.vwa.small_kd_min = 0.3;

cfg.mf12 = struct();
cfg.mf12.energy_keep = 0.995;
cfg.mf12.max_components = 256;
cfg.mf12.drop_zero_frequency = true;

cfg.validation = struct();
cfg.validation.enable_templates_only = true;
end
