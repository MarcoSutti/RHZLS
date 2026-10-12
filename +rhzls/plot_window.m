function [window, k_converged] = plot_window(ax, history)
%PLOT_WINDOW Show all recorded iterations, including the numerical plateaus.
% The first |grad| <= 1e-14 crossing is diagnostic only and does not crop data.
% Markers follow the complete viewport.
lines = findobj(ax, 'Type', 'line');
last = 0;
for h = reshape(lines,1,[])
    values = h.XData;
    if ~isempty(values), last = max(last, max(values)); end
end
k_converged = NaN;
index = find([history.gradnorm] <= 1e-14, 1, 'first');
if ~isempty(index)
    k_converged = history(index).iter;
end
window = [0, max(1,last)];
xlim(ax,window);
for h = reshape(lines,1,[])
    if strcmp(h.Marker,'none'), continue; end
    visible = find(h.XData >= window(1) & h.XData <= window(2));
    if ~isempty(visible)
        h.MarkerIndices = visible(unique(round(linspace(1,numel(visible),min(12,numel(visible))))));
    end
end
end
