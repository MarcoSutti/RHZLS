function [y, velocity] = oblique_curve(x, d, t)
%OBLIQUE_CURVE Column normalization and exact derivative on a product of spheres.
if nargin < 3, t = 1; end
if ~isscalar(t) || ~isreal(t) || ~isfinite(t) || ...
        ~isnumeric(x) || ~isnumeric(d) || ~ismatrix(x) || isempty(x) || ...
        ~isequal(size(x),size(d)) || ~isreal(x) || ~isreal(d) || ...
        any(~isfinite(x(:))) || any(~isfinite(d(:)))
    error('rhzls:InvalidRetraction','Oblique retraction requires finite real compatible inputs.');
end
z=x+t*d; scale=sqrt(sum(z.^2,1));
if any(~isfinite(scale)) || any(scale==0)
    error('rhzls:InvalidRetraction','Invalid column normalization.');
end
y=bsxfun(@rdivide,z,scale);
if nargout>1
    velocity=bsxfun(@rdivide,d-bsxfun(@times,y,sum(y.*d,1)),scale);
end
end
