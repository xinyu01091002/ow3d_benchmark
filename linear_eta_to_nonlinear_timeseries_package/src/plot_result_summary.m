function figs = plot_result_summary(result, output_dir, cfg)
%PLOT_RESULT_SUMMARY Create standard summary plots.

if nargin < 3 || isempty(cfg)
    cfg = default_config();
end

t = result.time(:);
Tp = result.meta.Tp;
[~, peak_idx] = max(abs(result.linear.eta(:)));
t_norm = (t - t(peak_idx)) / Tp;
figs = struct();

figs.timeseries = figure('Color', 'w', 'Position', [80, 80, 1200, 900], 'Visible', 'off');
tiledlayout(3, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(t_norm, result.linear.eta, 'k-', 'LineWidth', 1.2, 'DisplayName', 'linear'); hold on;
plot(t_norm, result.superharmonic.eta, 'r--', 'LineWidth', 1.2, 'DisplayName', 'superharmonic');
plot(t_norm, result.subharmonic.eta, 'b-.', 'LineWidth', 1.2, 'DisplayName', 'subharmonic');
plot(t_norm, result.nonlinear.eta, 'g-', 'LineWidth', 1.4, 'DisplayName', 'nonlinear');
ylabel('\eta (m)');
title('Surface Elevation Time Series');
xlim([-3, 3]);
grid on;
legend('Location', 'best');

nexttile;
plot(t_norm, result.linear.u_surface, 'k-', 'LineWidth', 1.2, 'DisplayName', 'linear'); hold on;
plot(t_norm, result.superharmonic.u_surface, 'r--', 'LineWidth', 1.2, 'DisplayName', 'superharmonic');
plot(t_norm, result.subharmonic.u_surface, 'b-.', 'LineWidth', 1.2, 'DisplayName', 'subharmonic');
plot(t_norm, result.nonlinear.u_surface, 'g-', 'LineWidth', 1.4, 'DisplayName', 'nonlinear');
ylabel('u_s (m/s)');
title('Surface Horizontal Velocity');
xlim([-3, 3]);
grid on;
legend('Location', 'best');

nexttile;
plot(t_norm, result.linear.w_surface, 'k-', 'LineWidth', 1.2, 'DisplayName', 'linear'); hold on;
plot(t_norm, result.superharmonic.w_surface, 'r--', 'LineWidth', 1.2, 'DisplayName', 'superharmonic');
plot(t_norm, result.subharmonic.w_surface, 'b-.', 'LineWidth', 1.2, 'DisplayName', 'subharmonic');
plot(t_norm, result.nonlinear.w_surface, 'g-', 'LineWidth', 1.4, 'DisplayName', 'nonlinear');
xlabel('(t - t_{peak}^{(1)}) / T_p');
ylabel('w_s (m/s)');
title('Surface Vertical Velocity');
xlim([-3, 3]);
grid on;
legend('Location', 'best');

figs.spectra = figure('Color', 'w', 'Position', [140, 140, 1200, 500], 'Visible', 'off');
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

dt = mean(diff(t));
omega = 2 * pi * [0:ceil(numel(t) / 2) - 1, -floor(numel(t) / 2):-1]' / (numel(t) * dt);

nexttile;
plot(fftshift(omega), fftshift(abs(fft(result.linear.eta)) / numel(t)), 'k-', 'LineWidth', 1.2); hold on;
plot(fftshift(omega), fftshift(abs(fft(result.superharmonic.eta)) / numel(t)), 'r--', 'LineWidth', 1.2);
plot(fftshift(omega), fftshift(abs(fft(result.subharmonic.eta)) / numel(t)), 'b-.', 'LineWidth', 1.2);
    plot(fftshift(omega), fftshift(abs(fft(result.nonlinear.eta)) / numel(t)), 'g-', 'LineWidth', 1.4);
xlabel('\omega (rad/s)');
ylabel('|FFT(\eta)|');
title('Elevation Spectra');
grid on;

nexttile;
closure_eta = result.nonlinear.eta - (result.linear.eta + result.superharmonic.eta + result.subharmonic.eta);
closure_u = result.nonlinear.u_surface - (result.linear.u_surface + result.superharmonic.u_surface + result.subharmonic.u_surface);
closure_w = result.nonlinear.w_surface - (result.linear.w_surface + result.superharmonic.w_surface + result.subharmonic.w_surface);
plot(t_norm, closure_eta, 'k-', 'LineWidth', 1.2, 'DisplayName', 'eta closure'); hold on;
plot(t_norm, closure_u, 'r--', 'LineWidth', 1.2, 'DisplayName', 'u closure');
plot(t_norm, closure_w, 'b-.', 'LineWidth', 1.2, 'DisplayName', 'w closure');
xlabel('(t - t_{peak}^{(1)}) / T_p');
ylabel('closure');
title('Assembly Closure Check');
xlim([-3, 3]);
grid on;
legend('Location', 'best');

if cfg.save_figures
    exportgraphics(figs.timeseries, fullfile(output_dir, 'summary_timeseries.png'), 'Resolution', 180);
    exportgraphics(figs.spectra, fullfile(output_dir, 'summary_spectra_and_closure.png'), 'Resolution', 180);
end
end
