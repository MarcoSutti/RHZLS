function handles = plot_history(name, info, reference, axes_pair, tolerance)
%PLOT_HISTORY Draw objective and gradient history with the shared style.
% AXES_PAIR contains the objective axes followed by the gradient axes.
if nargin < 5, tolerance = []; end
if isempty(info), handles = gobjects(0); return; end
style = rhzls.plot_style(name, numel(info));
ax_cost = axes_pair(1); ax_grad = axes_pair(2);
if isfield(reference, 'optimum')
    values = abs([info.cost] - reference.optimum);
    cost_scale = abs(reference.optimum);
    cost_label = '$|f(x_k)-f_\star|/|f_\star|$';
    if cost_scale == 0
        cost_scale = values(1);
        cost_label = '$|f(x_k)-f_\star|/|f(x_0)-f_\star|$';
    end
    if cost_scale == 0, cost_scale = 1; end
    values = values / cost_scale;
    values(values == 0) = eps;
    handles(1) = semilogy(ax_cost, [info.iter], values, style.line_args{:});
    ylabel(ax_cost, cost_label, 'Interpreter', 'latex');
else
    handles(1) = plot(ax_cost, [info.iter], [info.cost], style.line_args{:});
    ylabel(ax_cost, '$f(x_k)$', 'Interpreter', 'latex');
end
gradient_scale = info(1).gradnorm;
if gradient_scale == 0, gradient_scale = 1; end
handles(2) = semilogy(ax_grad, [info.iter], max([info.gradnorm]/gradient_scale, eps), style.line_args{:});
ylabel(ax_grad, '$\|g_k\|/\|g_0\|$', 'Interpreter', 'latex');
set(axes_pair, 'FontSize', style.font_size, ...
    'TickLabelInterpreter', 'latex', 'GridColor', style.grid_color);
for ax = axes_pair(:).'
    xlabel(ax, 'iteration $k$', 'Interpreter', 'latex');
    grid(ax, 'on');
    rhzls.plot_window(ax,info);
end
if ~isempty(tolerance) && tolerance > 0
    yline(ax_grad, tolerance/gradient_scale, '--', 'tolerance', ...
        'Color', style.reference_color, 'LineWidth', 2, 'Interpreter', 'latex');
end
end
