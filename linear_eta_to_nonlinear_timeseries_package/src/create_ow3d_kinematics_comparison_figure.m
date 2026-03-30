function metrics = create_ow3d_kinematics_comparison_figure(output_dir, ow3d, result)
%CREATE_OW3D_KINEMATICS_COMPARISON_FIGURE Compare eta/phi/u/w totals against OW3D probe series.

t = ow3d.time(:);
[~, peak_idx] = max(abs(ow3d.eta.first));
Tp = result.meta.Tp;
t_plot = (t - t(peak_idx)) / Tp;
case_text = build_case_text_local(ow3d, result);

metrics = struct();
metrics.eta = compute_metrics_local(ow3d.eta.total(:), result.nonlinear.eta(:));
metrics.phi = compute_metrics_local(ow3d.phi.total(:), result.nonlinear.phi_surface(:));
metrics.u = compute_metrics_local(ow3d.u.total(:), result.nonlinear.u_surface(:));
metrics.w = compute_metrics_local(ow3d.w.total(:), result.nonlinear.w_surface(:));

fig = figure('Color', 'w', 'Position', [100, 100, 1200, 900], 'Visible', 'off');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
sgtitle(sprintf('OW3D kinematics probe vs VWA package (%s)', case_text), 'Interpreter', 'none');
draw_panel_local(nexttile, t_plot, ow3d.eta.total(:), result.nonlinear.eta(:), ...
    sprintf('eta nonlinear: corr=%.4f rmse=%.3e', metrics.eta.corr, metrics.eta.rmse), '\eta (m)');
draw_panel_local(nexttile, t_plot, ow3d.phi.total(:), result.nonlinear.phi_surface(:), ...
    sprintf('phi nonlinear: corr=%.4f rmse=%.3e', metrics.phi.corr, metrics.phi.rmse), '\phi_s (m^2/s)');
draw_panel_local(nexttile, t_plot, ow3d.u.total(:), result.nonlinear.u_surface(:), ...
    sprintf('u nonlinear: corr=%.4f rmse=%.3e', metrics.u.corr, metrics.u.rmse), 'u_s (m/s)');
draw_panel_local(nexttile, t_plot, ow3d.w.total(:), result.nonlinear.w_surface(:), ...
    sprintf('w nonlinear: corr=%.4f rmse=%.3e', metrics.w.corr, metrics.w.rmse), 'w_s (m/s)');
xlabel('(t - t_{peak}^{(1)}) / T_p');
exportgraphics(fig, fullfile(output_dir, 'ow3d_kinematics_probe_comparison.png'), 'Resolution', 180);
close(fig);
save(fullfile(output_dir, 'ow3d_kinematics_probe_metrics.mat'), 'metrics');
end

function draw_panel_local(ax, t_plot, ref, cand, title_text, ylabel_text)
plot(ax, t_plot, ref, 'k-', 'LineWidth', 1.3, 'DisplayName', 'OW3D'); hold(ax, 'on');
plot(ax, t_plot, cand, 'r--', 'LineWidth', 1.3, 'DisplayName', 'VWA package');
xlim(ax, [-3, 3]);
grid(ax, 'on');
title(ax, title_text);
ylabel(ax, ylabel_text);
legend(ax, 'Location', 'best');
end

function m = compute_metrics_local(ref, cand)
ref = ref(:); cand = cand(:);
m = struct();
m.rmse = sqrt(mean((ref-cand).^2));
if all(abs(ref-mean(ref)) < eps) || all(abs(cand-mean(cand)) < eps)
    m.corr = NaN;
else
    c = corrcoef(ref, cand);
    m.corr = c(1,2);
end
m.peak_ratio = max(abs(cand)) / max(abs(ref)+eps);
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
    if isfield(ow3d, 'probe_x')
        parts(end+1) = sprintf('x=%.3f m', ow3d.probe_x);
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
