function [y, velocity] = generalized_grassmann_curve(x, d, t, L)
%GENERALIZED_GRASSMANN_CURVE Polar retraction for X'*B*X=I, B=L*L'.
% L is the lower Cholesky factor of B. Differentiate in whitened coordinates.
if nargout > 1
    [q, dq] = rhzls.retractions.stiefel_curve(L'*x, L'*d, t);
    y = L'\q;
    velocity = L'\dq;
else
    y = L'\rhzls.retractions.stiefel_curve(L'*x, L'*d, t);
end
end
