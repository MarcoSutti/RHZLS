function position = figure_position()
%FIGURE_POSITION Centered window occupying approximately half the screen area.
% Both dimensions use the same factor, preserving the screen aspect ratio.
previous_units = get(groot, 'Units');
restore_units = onCleanup(@() set(groot, 'Units', previous_units)); %#ok<NASGU>
set(groot, 'Units', 'pixels');
screen = get(groot, 'ScreenSize');
size_pixels = round(screen(3:4) / sqrt(2));
origin = screen(1:2) + round((screen(3:4) - size_pixels) / 2);
position = [origin, size_pixels];
end
