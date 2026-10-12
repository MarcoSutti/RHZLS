function [y, velocity] = stiefel_curve(x, d, t)
% Compatibility entry point. Implementation: +rhzls/+retractions/stiefel_curve.m
if nargin < 3, t = 1; end
if nargout > 1
    [y, velocity] = rhzls.retractions.stiefel_curve(x, d, t);
else
    y = rhzls.retractions.stiefel_curve(x, d, t);
end
end