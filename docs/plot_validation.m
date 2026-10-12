% Reproduce and show validation plots from saved results, without optimization.
% Interactive MATLAB displays the figures; batch runs can set
% set(groot, 'defaultFigureVisible', 'off') to export without opening windows.
folder = fileparts(mfilename('fullpath'));
show_figures = strcmp(char(get(groot, 'defaultFigureVisible')), 'on');
data = load(fullfile(folder, 'demo_results.mat'), 'info', 'details');
names = {'stiefel', 'sphere'};
labels = {'Brockett / Stiefel', 'Rayleigh / sphere'};
if isfield(data.details.stiefel, 'comparison')
    rhzls.plot_style('stiefel', numel(data.info.stiefel));
    for k = 1:numel(names)
        comparison = data.details.(names{k}).comparison;
        fig = rhzls.plot_comparison(comparison);
        rhzls.export_pdf(fig, fullfile(folder, ['comparison_' names{k} '.pdf']));
        if ~show_figures, close(fig); end
    end
    if exist(fullfile(folder,'accuracy_examples.mat'),'file')
        plot_accuracy_examples;
    end
    if exist(fullfile(folder,'manopt_candidate_results.mat'),'file')
        plot_manopt_examples;
    end
    return;
end
rhzls.plot_style('stiefel', numel(data.info.stiefel));
fig = figure('Color', 'w', 'WindowStyle', 'normal', 'Units', 'pixels', 'Position', rhzls.figure_position());
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:numel(names)
    history = data.info.(names{k});
    axes_pair = [nexttile(layout, k), nexttile(layout, k+numel(names))];
    rhzls.plot_history(names{k}, history, data.details.(names{k}).reference, ...
        axes_pair, data.details.(names{k}).options.tolgradnorm);
    title(axes_pair(1), labels{k}, 'Interpreter', 'latex');
end
title(layout, 'RHZLS: one line search, two spectral examples (seed 10)');
rhzls.export_pdf(fig, fullfile(folder, 'convergence.pdf'));
if ~show_figures, close(fig); end
