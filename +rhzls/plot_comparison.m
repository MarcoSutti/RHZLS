function [fig, handles] = plot_comparison(result, ax)
%PLOT_COMPARISON Normalized convergence and line-search cost evaluations.
% Solid: relative objective gap; dashed: normalized gradient; dotted: relative solution error.
% Numerical zeros are displayed at eps on the log axis; histories stay raw.
% The lower panel shows exact cost callback counts per completed line search.
% A supplied axes retains the convergence-only interface for embedded plots.
project = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project, 'options_plot.m'));
two_panels = nargin < 2 || isempty(ax);
methods = {'manopt', 'rhzls'};
method_names = {'Manopt (Armijo)', 'RHZLS'};
if isfield(result, 'rhzls_label'), method_names{2} = result.rhzls_label; end
if ~isfield(result, 'rhzls_label') && isfield(result.rhzls, 'details') ...
        && isfield(result.rhzls.details.options, 'rhzls') ...
        && isfield(result.rhzls.details.options.rhzls, 'acceptance') ...
        && strcmp(result.rhzls.details.options.rhzls.acceptance, 'approximate')
    method_names{2} = 'RHZLS (approx. Wolfe)';
end
colors = {Olive, Candy}; markers = {'o', '^'};
if isfield(result, 'methods')
    methods = result.methods;
    method_names = result.method_labels;
    colors = {Olive, gray3, Candy}; markers = {'o', 's', '^'};
end
if two_panels
    for m = 1:numel(methods)
        if ~isfield(result.(methods{m}).info, 'nf')
            error('rhzls:EvaluationCountsUnavailable', ...
                ['Exact line-search cost counts are absent for %s. ' ...
                 'Run run_all_experiments(false) to regenerate the results.'], method_names{m});
        end
    end
    fig = figure('Name', ['RHZLS vs Manopt / ' result.name], ...
        'NumberTitle', 'off', 'Color', 'w', 'WindowStyle', 'normal', ...
        'Units', 'pixels', 'Position', rhzls.figure_position());
    % Leave room around the external legends, including long LaTeX labels.
    layout = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'loose');
    % Reserve an extra right margin for LaTeX glyphs in the PDF export.
    layout.OuterPosition = [0, 0, 0.97, 1];
    ax = nexttile(layout, 1, [2, 1]);
    evaluations_ax = nexttile(layout, 3);
else
    fig = ancestor(ax, 'figure');
end
cla(ax); hold(ax, 'on');
handles = gobjects(0); labels = {};
initial = result.(methods{1}).info(1);
cost_scale = abs(result.reference.optimum);
cost_label = '$|f_k-f_\star|/|f_\star|$';
if cost_scale == 0
    cost_scale = abs(initial.cost - result.reference.optimum);
    cost_label = '$|f_k-f_\star|/|f_0-f_\star|$';
end
if cost_scale == 0, cost_scale = 1; end
gradient_scale = initial.gradnorm;
if gradient_scale == 0, gradient_scale = 1; end
fields = {'objective_error', 'gradient_error', 'solution_error'};
metric_labels = {cost_label, ...
    '$\|g_k\|/\|g_0\|$', result.reference.solution_label};
line_styles = {'-', '--', ':'};
for m = 1:numel(methods)
    history = result.(methods{m}).info;
    indices = unique(round(linspace(1, numel(history), min(12, numel(history)))));
    for j = 1:3
        if j == 1
            % Normalize saved absolute gaps without changing raw histories.
            values = abs([history.cost] - result.reference.optimum);
            if isfield(result.reference, 'direct_objective_error') && result.reference.direct_objective_error
                values = [history.objective_error];
            end
            values = values / cost_scale;
        elseif j == 2
            values = [history.gradnorm] / gradient_scale;
        else
            values = [history.(fields{j})];
        end
        if all(isnan(values)), continue; end
        values(values == 0) = eps;
        h = semilogy(ax, [history.iter], values, ...
            'Color', colors{m}, 'LineStyle', line_styles{j}, 'LineWidth', 2, ...
            'Marker', markers{m}, 'MarkerIndices', indices, 'MarkerSize', 10, ...
            'MarkerFaceColor', colors{m}, 'MarkerEdgeColor', 'k');
        handles(end+1) = h; %#ok<AGROW>
        labels{end+1} = [metric_labels{j} ', ' method_names{m}]; %#ok<AGROW>
    end
end
yline(ax, sqrt(eps), '--', 'Color', gray1, 'LineWidth', 2, 'HandleVisibility', 'off');
yline(ax, eps, '--', 'Color', gray1, 'LineWidth', 2, 'HandleVisibility', 'off');
set(ax, 'YScale', 'log', 'FontSize', get(groot, 'defaultAxesFontSize'), ...
    'TickLabelInterpreter', 'latex', 'GridColor', gray1, 'Box', 'on');
grid(ax, 'on');
xlabel(ax, 'iteration $k$ of steepest descent', 'Interpreter', 'latex');
ylabel(ax, 'normalized convergence errors', 'Interpreter', 'latex');
switch result.name
    case 'stiefel', plot_title = 'Brockett / Stiefel';
    case 'sphere', plot_title = 'Rayleigh / sphere';
    otherwise, plot_title = result.name;
end
if isfield(result, 'display_title'), plot_title = result.display_title; end
title(ax, plot_title, 'Interpreter', 'latex');
legend_size = 16;
if numel(methods) > 2, legend_size = 12; end
legend(ax, handles, labels, 'Interpreter', 'latex', ...
    'FontSize', legend_size, 'Location', 'northeastoutside', 'Box', 'off');
hold(ax, 'off');
window_history = result.rhzls.info;
if isfield(result, 'window_history'), window_history = result.window_history; end
[window, ~] = rhzls.plot_window(ax, window_history);
if two_panels
    % Keep both panels legible in the fixed half-screen figure window.
    set(ax, 'FontSize', 16);
    set([ax.Title, ax.YLabel], 'FontSize', 17);
    xlabel(ax, '');
    set(legend(ax), 'FontSize', 12);
    set(handles, 'MarkerSize', 8);
    rhzls.plot_evaluations(result, evaluations_ax, methods, method_names, colors, window, gray1);
    xlim(evaluations_ax, window);
    linkaxes([ax, evaluations_ax], 'x');
    % Keep the public CurrentAxes convention and the six convergence handles.
    set(fig, 'CurrentAxes', ax);
else
    % Compact typography for panels sharing the common figure window.
    set(ax, 'FontSize', 12);
    set([ax.Title, ax.XLabel, ax.YLabel], 'FontSize', 12);
    set(legend(ax), 'FontSize', 9);
    set(handles, 'MarkerSize', 6);
end
end
