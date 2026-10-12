function style = plot_style(name, point_count)
%PLOT_STYLE Load the user's options_plot.m and select a curve style.
% Colors and LaTeX/font defaults come from RHZLS/options_plot.m.
% About 12 filled markers include both the initial and final iteration.
project = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(project, 'options_plot.m'));
switch lower(char(name))
    case 'stiefel'
        color = Orchid; marker = 'o';
    case 'sphere'
        color = Cerulean; marker = 's';
    case 'fixedrank'
        color = Candy; marker = 'd';
    otherwise
        color = blue; marker = 'o';
end
indices = unique(round(linspace(1, point_count, min(12, point_count))));
style = struct('font_size', get(groot, 'defaultAxesFontSize'), ...
    'grid_color', gray1, 'reference_color', gray2, ...
    'line_args', {{'Color', color, 'LineWidth', 2, 'Marker', marker, ...
        'MarkerSize', 9, 'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', ...
        'MarkerIndices', indices}});
end
