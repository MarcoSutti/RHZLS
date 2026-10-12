function [y, velocity] = stiefel_curve(x, d, t)
%STIEFEL_CURVE Polar retraction and its exact curve velocity.
%   Write z=x+t*d=y*H with H positive definite. Differentiate H^2=z'*z:
%   H*Hdot+Hdot*H=d'*z+z'*d, then velocity=(d-y*Hdot)/H.
%   The Sylvester equation is solved in H's eigenbasis, so no additional
%   toolbox is required. This formula also accounts for roundoff in the
%   tangent constraint x'*d+d'*x=0.
if nargin < 3
    t = 1;
end
if ~isscalar(t) || ~isreal(t) || ~isfinite(t) || ...
        ~isnumeric(x) || ~isnumeric(d) || ~ismatrix(x) || ...
        size(x, 2) == 0 || size(x, 1) < size(x, 2) || ...
        ~isequal(size(x), size(d)) || ~isreal(x) || ~isreal(d) || ...
        any(~isfinite(x(:))) || any(~isfinite(d(:)))
    error('rhzls:InvalidRetraction', 'Stiefel retraction requires finite real compatible matrices.');
end
z = x + t*d;
if any(~isfinite(z(:)))
    error('rhzls:InvalidRetraction', 'Nonfinite polar retraction trial.');
end
[u, singularvalues, q] = svd(z, 'econ');
s = diag(singularvalues);
if any(~isfinite(s)) || any(s <= 0)
    error('rhzls:InvalidRetraction', 'The polar retraction requires a full-column-rank trial.');
end
% The right singular vectors diagonalize H, without forming z'*z and
% squaring its condition number.
y = u*q.';
if nargout > 1
    rhs = d.'*z + z.'*d;
    hdot_basis = (q.'*rhs*q)./bsxfun(@plus, s, s.');
    hdot = q*hdot_basis*q.';
    velocity = bsxfun(@rdivide, (d-y*hdot)*q, s.')*q.';
    if any(~isfinite(velocity(:)))
        error('rhzls:InvalidRetraction', 'Nonfinite polar retraction velocity.');
    end
end
end
