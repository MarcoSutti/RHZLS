function export_pdf(fig, filename)
%EXPORT_PDF Vector PDF with the same page size for every plot.
position = rhzls.figure_position();
page_inches = position(3:4) / get(groot, 'ScreenPixelsPerInch');
set(fig, 'PaperUnits', 'inches', 'PaperSize', page_inches, ...
    'PaperPosition', [0 0 page_inches], 'PaperPositionMode', 'manual');
drawnow;
% The vector PDF output preserves transparent fills and LaTeX tick labels.
print(fig, filename, '-dpdf', '-vector');
end
