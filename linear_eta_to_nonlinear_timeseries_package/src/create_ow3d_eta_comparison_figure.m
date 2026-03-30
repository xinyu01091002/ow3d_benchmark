function metrics = create_ow3d_eta_comparison_figure(output_dir, ow3d, result)
%CREATE_OW3D_ETA_COMPARISON_FIGURE Compare package eta outputs with OW3D harmonics.

t = ow3d.time(:);
[~, peak_idx] = max(abs(ow3d.eta.first));
t_center = t(peak_idx);
Tp = get_tp_local(ow3d, result);
t_plot = (t - t_center) / Tp;
case_text = build_case_text_local(ow3d, result);

ow_eta_linear = ow3d.eta.first(:);
ow_eta_super = ow3d.eta.second_super(:) + ow3d.eta.third(:);
ow_eta_sub = ow3d.eta.second_sub(:);
ow_eta_total = ow3d.eta.total(:);

pkg_eta_linear = result.linear.eta(:);
pkg_eta_super = result.superharmonic.eta(:);
pkg_eta_sub = result.subharmonic.eta(:);
pkg_eta_total = result.nonlinear.eta(:);

metrics = struct();
metrics.linear = compute_pair_metrics_local(ow_eta_linear, pkg_eta_linear);
metrics.super = compute_pair_metrics_local(ow_eta_super, pkg_eta_super);
metrics.sub = compute_pair_metrics_local(ow_eta_sub, pkg_eta_sub);
metrics.nonlinear = compute_pair_metrics_local(ow_eta_total, pkg_eta_total);

fig = figure('Color', 'w', 'Position', [100, 100, 1200, 900], 'Visible', 'off');
tiledlayout(4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(sprintf('OW3D vs VWA package eta comparison (%s)', case_text), 'Interpreter', 'none');

draw_panel_local(nexttile, t_plot, ow_eta_linear, pkg_eta_linear, ...
    sprintf('Linear eta: corr=%.4f, rmse=%.3e', metrics.linear.corr, metrics.linear.rmse));
draw_panel_local(nexttile, t_plot, ow_eta_super, pkg_eta_super, ...
    sprintf('Superharmonic eta: corr=%.4f, rmse=%.3e', metrics.super.corr, metrics.super.rmse));
draw_panel_local(nexttile, t_plot, ow_eta_sub, pkg_eta_sub, ...
    sprintf('Subharmonic eta: corr=%.4f, rmse=%.3e', metrics.sub.corr, metrics.sub.rmse));
draw_panel_local(nexttile, t_plot, ow_eta_total, pkg_eta_total, ...
    sprintf('Nonlinear eta: corr=%.4f, rmse=%.3e', metrics.nonlinear.corr, metrics.nonlinear.rmse));
xlabel('(t - t_{peak}^{(1)}) / T_p');

exportgraphics(fig, fullfile(output_dir, 'ow3d_eta_comparison.png'), 'Resolution', 180);
close(fig);

save(fullfile(output_dir, 'ow3d_eta_comparison_metrics.mat'), 'metrics');
end

function draw_panel_local(ax, t_plot, ow_vals, pkg_vals, title_text)
plot(ax, t_plot, ow_vals, 'k-', 'LineWidth', 1.3, 'DisplayName', 'OW3D'); hold(ax, 'on');
plot(ax, t_plot, pkg_vals, 'r--', 'LineWidth', 1.3, 'DisplayName', 'VWA package');
xlim(ax, [-3, 3]);
grid(ax, 'on');
ylabel(ax, '\eta (m)');
title(ax, title_text);
legend(ax, 'Location', 'best');
end

function metrics = compute_pair_metrics_local(ref_vals, test_vals)
ref_vals = ref_vals(:);
test_vals = test_vals(:);
metrics = struct();
metrics.rmse = sqrt(mean((ref_vals - test_vals).^2));
if all(abs(ref_vals - mean(ref_vals)) < eps) || all(abs(test_vals - mean(test_vals)) < eps)
    metrics.corr = NaN;
else
    c = corrcoef(ref_vals, test_vals);
    metrics.corr = c(1, 2);
end
metrics.peak_ratio = max(abs(test_vals)) / max(abs(ref_vals) + eps);
end

function Tp = get_tp_local(~, result)
Tp = result.meta.Tp;
end

function txt = build_case_text_local(ow3d, result)
parts = strings(0, 1);
if isfield(ow3d, 'case_cfg')
    cfg = ow3d.case_cfg;
    if isfield(cfg, 'kd')
        parts(end+1) = sprintf('kd=%.3g', cfg.kd);
    end
    if isfield(cfg, 'Alpha')
        parts(end+1) = sprintf('Alpha=%.3g', cfg.Alpha);
    end
    if isfield(cfg, 'Akp')
        parts(end+1) = sprintf('Akp=%.3g', cfg.Akp);
    end
end
if isempty(parts) && isfield(result, 'meta')
    parts(end+1) = sprintf('h=%.3f m', result.meta.h);
end
txt = strjoin(cellstr(parts), ', ');
if strlength(txt) == 0
    txt = 'case metadata unavailable';
end
end
