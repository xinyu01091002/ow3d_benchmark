function plot_bingchen_test_case_overview()
%PLOT_BINGCHEN_TEST_CASE_OVERVIEW Create quick-look overview plots for Bingchen test cases.

this_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(fileparts(this_dir));
input_dir = fullfile(root_dir, 'examples', 'bingchen_test_cases');
output_dir = fullfile(root_dir, 'output_example', 'bingchen_test_cases');
case_files = dir(fullfile(input_dir, 'case_*.mat'));
n_cases = numel(case_files);
if n_cases == 0
    error('No Bingchen cases found in %s', input_dir);
end

case_order = zeros(n_cases, 1);
for idx = 1:n_cases
    tmp = load(fullfile(case_files(idx).folder, case_files(idx).name), 'case_meta');
    case_order(idx) = tmp.case_meta.Tp;
end
[~, order_idx] = sort(case_order);
case_files = case_files(order_idx);

clr_input = [0.08, 0.32, 0.58];
clr_linear = [0.12, 0.12, 0.12];
clr_nonlinear = [0.00, 0.50, 0.30];
font_name = 'Times New Roman';

fig_input = figure('Color', 'w', 'Position', [60, 60, 1400, 900], 'Visible', 'off');
tile_input = tiledlayout(n_cases, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(tile_input, 'Bingchen Liu Test Cases: Input Linear \eta(t) Centered At Main Peak', ...
    'FontName', font_name, 'FontSize', 18, 'FontWeight', 'bold');

fig_output = figure('Color', 'w', 'Position', [80, 80, 1400, 1200], 'Visible', 'off');
tile_output = tiledlayout(n_cases, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
sgtitle(tile_output, 'Bingchen Liu Test Cases: Linear Vs Nonlinear Package Response', ...
    'FontName', font_name, 'FontSize', 18, 'FontWeight', 'bold');

for idx = 1:n_cases
    input_path = fullfile(case_files(idx).folder, case_files(idx).name);
    input_data = load(input_path);
    [~, case_name] = fileparts(case_files(idx).name);
    result_path = fullfile(output_dir, case_name, 'result_bundle.mat');
    if ~isfile(result_path)
        error('Missing result bundle: %s', result_path);
    end
    result_data = load(result_path);
    result = result_data.result;

    t = input_data.t(:);
    Tp = input_data.case_meta.Tp;
    [~, input_peak_idx] = max(abs(input_data.eta_linear(:)));
    t_peak_centered = (t - t(input_peak_idx)) / Tp;
    case_label = sprintf('%s  |  T_p=%.3g s, A=%.3g m, k_p=%.3g 1/m, h=%.3g m', ...
        case_name, input_data.case_meta.Tp, input_data.case_meta.A, input_data.case_meta.kp, input_data.case_meta.h);

    figure(fig_input);
    ax = nexttile;
    plot(ax, t_peak_centered, input_data.eta_linear(:), '-', 'Color', clr_input, 'LineWidth', 2.1);
    hold(ax, 'on');
    xline(ax, 0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.1);
    yline(ax, 0, '-', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);
    grid(ax, 'on');
    xlim(ax, [-3, 3]);
    ylabel(ax, '\eta_{linear} (m)');
    title(ax, case_label, 'Interpreter', 'none');
    set(ax, 'FontName', font_name, 'FontSize', 12, 'LineWidth', 1.0, 'Box', 'on');
    if idx == n_cases
        xlabel(ax, '(t - t_{peak}^{(1)}) / T_p');
    end

    [~, peak_idx] = max(abs(result.linear.eta(:)));
    t_peak_norm = result.time(:) / result.meta.Tp - result.time(peak_idx) / result.meta.Tp;

    figure(fig_output);
    ax1 = nexttile;
    plot(ax1, t_peak_norm, result.linear.eta(:), '-', 'Color', clr_linear, 'LineWidth', 1.7, 'DisplayName', 'linear'); hold(ax1, 'on');
    plot(ax1, t_peak_norm, result.nonlinear.eta(:), '-', 'Color', clr_nonlinear, 'LineWidth', 2.1, 'DisplayName', 'nonlinear');
    style_panel_local(ax1, [-3, 3], font_name);
    ylabel(ax1, '\eta (m)');
    title(ax1, sprintf('%s | \\eta', case_label), 'Interpreter', 'none');
    legend(ax1, 'Location', 'northwest', 'Box', 'off');

    ax2 = nexttile;
    plot(ax2, t_peak_norm, result.linear.u_surface(:), '-', 'Color', clr_linear, 'LineWidth', 1.7, 'DisplayName', 'linear'); hold(ax2, 'on');
    plot(ax2, t_peak_norm, result.nonlinear.u_surface(:), '-', 'Color', clr_nonlinear, 'LineWidth', 2.1, 'DisplayName', 'nonlinear');
    style_panel_local(ax2, [-3, 3], font_name);
    ylabel(ax2, 'u_s (m/s)');
    title(ax2, sprintf('%s | u_s', case_label), 'Interpreter', 'none');
    legend(ax2, 'Location', 'northwest', 'Box', 'off');

    ax3 = nexttile;
    plot(ax3, t_peak_norm, result.linear.w_surface(:), '-', 'Color', clr_linear, 'LineWidth', 1.7, 'DisplayName', 'linear'); hold(ax3, 'on');
    plot(ax3, t_peak_norm, result.nonlinear.w_surface(:), '-', 'Color', clr_nonlinear, 'LineWidth', 2.1, 'DisplayName', 'nonlinear');
    style_panel_local(ax3, [-3, 3], font_name);
    ylabel(ax3, 'w_s (m/s)');
    title(ax3, sprintf('%s | w_s', case_label), 'Interpreter', 'none');
    legend(ax3, 'Location', 'northwest', 'Box', 'off');
    if idx == n_cases
        xlabel(ax1, '(t - t_{peak}^{(1)}) / T_p');
        xlabel(ax2, '(t - t_{peak}^{(1)}) / T_p');
        xlabel(ax3, '(t - t_{peak}^{(1)}) / T_p');
    end
end

exportgraphics(fig_input, fullfile(output_dir, 'bingchen_input_eta_overview.png'), 'Resolution', 180);
exportgraphics(fig_output, fullfile(output_dir, 'bingchen_nonlinear_overview.png'), 'Resolution', 180);
close(fig_input);
close(fig_output);

fprintf('Bingchen overview plots written to:\n%s\n', output_dir);
end

function style_panel_local(ax, x_limits, font_name)
xline(ax, 0, ':', 'Color', [0.35, 0.35, 0.35], 'LineWidth', 1.1);
yline(ax, 0, '-', 'Color', [0.88, 0.88, 0.88], 'LineWidth', 0.8);
xlim(ax, x_limits);
grid(ax, 'on');
set(ax, 'FontName', font_name, 'FontSize', 11, 'LineWidth', 1.0, 'Box', 'on');
end
