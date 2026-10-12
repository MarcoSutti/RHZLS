function [y, velocity] = fixedrank_curve(x, d, t)
% Compatibility entry point. Implementation: +rhzls/+retractions/fixedrank_curve.m
if nargin < 3, t = 1; end
if nargout > 1
    [y, velocity] = rhzls.retractions.fixedrank_curve(x, d, t);
else
    y = rhzls.retractions.fixedrank_curve(x, d, t);
end
end