function handles = plot_evaluations(result, ax, methods, method_names, colors, window, grid_color)
%PLOT_EVALUATIONS Transparent filled steps of cost callbacks per line search.
% Each count n_f(k) occupies the interval (k-1,k]. Iteration zero is excluded.
% Counts end at each method's final iteration, without extending a last value.
% Methods, names and colors are supplied by plot_comparison for consistency.
cla(ax); hold(ax, 'on');
if nargin < 6, window = [-Inf, Inf]; end
if nargin < 7, grid_color = get(ax, 'GridColor'); end
handles = gobjects(0);
largest_count = 0;
for m = 1:numel(methods)
    history = result.(methods{m}).info;
    if ~isfield(history, 'nf')
        error('rhzls:EvaluationCountsUnavailable', ...
            ['Exact line-search cost counts are absent for %s. ' ...
             'Run run_all_experiments(false) to regenerate the results.'], method_names{m});
    end
    completed = [history.iter] > 0;
    iterations = [history(completed).iter];
    counts = [history(completed).nf];
    if isempty(iterations), continue; end
    if numel(counts) ~= numel(iterations) || any(~isfinite(counts)) ...
            || any(counts < 0) || any(counts ~= round(counts))
        error('rhzls:InvalidEvaluationCounts', ...
            'Line-search cost counts for %s must be nonnegative integers.', method_names{m});
    end
    if any(diff(iterations) ~= 1)
        error('rhzls:InvalidEvaluationHistory', ...
            'Completed iterations for %s must be consecutive.', method_names{m});
    end
    [step_x, step_y] = stairs([iterations(1)-1, iterations], [counts, counts(end)]);
    h = patch(ax, [step_x(1); step_x(:); step_x(end)], ...
        [0; step_y(:); 0], colors{m}, 'FaceAlpha', 0.4, ...
        'EdgeColor', 'k', 'LineWidth', 0.5, ...
        'DisplayName', method_names{m}, 'Tag', 'rhzlsFunctionEvaluations');
    handles(end+1) = h; %#ok<AGROW>
    visible = iterations > window(1) & iterations <= window(2);
    if any(visible), largest_count = max(largest_count, max(counts(visible))); end
end
set(ax, 'YScale', 'linear', 'FontSize', 16, 'TickLabelInterpreter', 'latex', ...
    'GridColor', grid_color, 'Box', 'on', 'Layer', 'top');
grid(ax, 'on');
upper = max(1, ceil(1.15*largest_count));
ylim(ax, [0, upper]);
yticks(ax, 0:max(1,ceil(upper/5)):upper);
xlabel(ax, 'iteration $k$ of steepest descent', 'Interpreter', 'latex', 'FontSize', 17);
ylabel(ax, '$n_f(k)$', 'Interpreter', 'latex', 'FontSize', 17);
if ~isempty(handles)
    legend(ax, handles, 'Interpreter', 'latex', 'FontSize', 12, ...
        'Location', 'northeastoutside', 'Box', 'off');
end
hold(ax, 'off');
end
