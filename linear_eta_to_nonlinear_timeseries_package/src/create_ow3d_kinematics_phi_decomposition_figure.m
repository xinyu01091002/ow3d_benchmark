function metrics = create_ow3d_kinematics_phi_decomposition_figure(output_dir, ow3d, result)
%CREATE_OW3D_KINEMATICS_PHI_DECOMPOSITION_FIGURE Compare phi decomposition against OW3D probe series.

t = ow3d.time(:);
[~, peak_idx] = max(abs(ow3d.eta.first));
Tp = result.meta.Tp;
t_plot = (t - t(peak_idx)) / Tp;
case_text = build_case_text_local(ow3d, result);

ref.linear = ow3d.phi.first(:);
ref.super = ow3d.phi.second_super(:) + ow3d.phi.third(:);
ref.sub = ow3d.phi.second_sub(:);
ref.nonlinear = ow3d.phi.total(:);

mod.linear = result.linear.phi_surface(:);
mod.super = result.superharmonic.phi_surface(:);
mod.sub = result.subharmonic.phi_surface(:);
mod.nonlinear = result.nonlinear.phi_surface(:);

labels = {'linear', 'super', 'sub', 'nonlinear'};
titles = struct( ...
    'linear', '\phi^{(1)}', ...
    'super', '\phi^{(2+)}+\phi^{(3)}', ...
    'sub', '\phi^{(2-)}', ...
    'nonlinear', '\phi nonlinear');

metrics = struct();
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 900], 'Visible', 'off');
tiledlayout(4,1,'Padding','compact','TileSpacing','compact');
sgtitle(sprintf('OW3D phi decomposition vs VWA package (%s)', case_text), 'Interpreter', 'none');
for i = 1:numel(labels)
    key = labels{i};
    metrics.(key) = compute_metrics_local(ref.(key), mod.(key), t);
    draw_panel_local(nexttile, t_plot, ref.(key), mod.(key), ...
        sprintf('%s: corr=%.4f rmse=%.3e lag=%d', titles.(key), metrics.(key).corr, metrics.(key).rmse, metrics.(key).peak_lag_samples), ...
        '\phi_s (m^2/s)');
end
xlabel('(t - t_{peak}^{(1)}) / T_p');
exportgraphics(fig, fullfile(output_dir, 'ow3d_kinematics_phi_decomposition.png'), 'Resolution', 180);
close(fig);

save(fullfile(output_dir, 'ow3d_kinematics_phi_decomposition_metrics.mat'), 'metrics');
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

function m = compute_metrics_local(ref, cand, t)
ref = ref(:);
cand = cand(:);
t = t(:);

m = struct();
m.rmse = sqrt(mean((ref - cand).^2));
if all(abs(ref - mean(ref)) < eps) || all(abs(cand - mean(cand)) < eps)
    m.corr = NaN;
else
    c = corrcoef(ref, cand);
    m.corr = c(1, 2);
end
m.peak_ratio = max(abs(cand)) / max(abs(ref) + eps);
[~, ref_idx] = max(abs(ref));
[~, cand_idx] = max(abs(cand));
m.reference_peak_time = t(ref_idx);
m.model_peak_time = t(cand_idx);
m.peak_lag_samples = cand_idx - ref_idx;
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
