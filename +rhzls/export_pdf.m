function export_pdf(fig, filename)
%EXPORT_PDF Vector PDF with the same page size for every plot.
position = rhzls.figure_position();
page_inches = position(3:4) / get(groot, 'ScreenPixelsPerInch');
set(fig, 'PaperUnits', 'inches', 'PaperSize', page_inches, ...
    'PaperPosition', [0 0 page_inches], 'PaperPositionMode', 'manual');
drawnow;
% A toolbar visible during interactive plotting can otherwise enter the PDF.
toolbars = findall(fig, 'Type', 'axestoolbar');
visibility = cell(size(toolbars));
visibility_modes = cell(size(toolbars));
for k = 1:numel(toolbars)
    visibility{k} = toolbars(k).Visible;
    if isprop(toolbars(k), 'VisibleMode')
        visibility_modes{k} = get(toolbars(k), 'VisibleMode');
    end
end
restore_toolbars = onCleanup(@() restore_toolbar_visibility( ...
    toolbars, visibility, visibility_modes)); %#ok<NASGU>
set(toolbars, 'Visible', 'off');
drawnow;
% The vector PDF output preserves transparent fills and LaTeX tick labels.
print(fig, filename, '-dpdf', '-vector');
end

function restore_toolbar_visibility(toolbars, visibility, visibility_modes)
for k = 1:numel(toolbars)
    if isvalid(toolbars(k))
        toolbars(k).Visible = visibility{k};
        if ~isempty(visibility_modes{k})
            set(toolbars(k), 'VisibleMode', visibility_modes{k});
        end
    end
end
end
