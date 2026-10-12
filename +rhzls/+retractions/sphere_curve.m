function [y, velocity] = sphere_curve(x, d, t)
%SPHERE_CURVE Normalized retraction and its exact derivative with respect to t.
%   y = (x+t*d)/norm(x+t*d), velocity = d/dt y.
if nargin < 3
    t = 1;
end
if ~isscalar(t) || ~isreal(t) || ~isfinite(t) || ...
        ~isnumeric(x) || ~isnumeric(d) || ~isequal(size(x), size(d)) || ...
        ~isreal(x) || ~isreal(d) || any(~isfinite(x(:))) || any(~isfinite(d(:)))
    error('rhzls:InvalidRetraction', 'Sphere retraction requires finite real compatible inputs.');
end
z = x + t*d;
scale = norm(z(:));
if ~isfinite(scale) || scale == 0
    error('rhzls:InvalidRetraction', 'The normalized sphere retraction is undefined at this trial.');
end
y = z/scale;
if nargout > 1
    velocity = (d - y*(y(:).'*d(:)))/scale;
end
end
