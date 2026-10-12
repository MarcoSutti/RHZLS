function bound = fixedrank_maxstep(x, d)
%FIXEDRANK_MAXSTEP Bound the positive smooth interval of the AO retraction.
%   S+t*M=S*(I+t*B), B=S\M, is singular at a positive t precisely when
%   B has a negative real eigenvalue lambda and t=-1/lambda. The bound
%   leaves a 10 percent margin before the first such pole, or is Inf when
%   there is no positive pole. A norm bound would be unnecessarily small
%   for nonnormal B. Nearly real eigenvalues are conservatively treated
%   as real; fixedrank_curve also guards numerical conditioning at every
%   actual trial. If eig fails, use the guaranteed local norm bound.
scaled = x.S\d.M;
magnitude = norm(scaled, 2);
if ~isfinite(magnitude)
    error('rhzls:InvalidRetraction', 'Cannot determine a finite orthographic retraction domain.');
end
if magnitude == 0
    bound = Inf;
    return;
end
try
    eigenvalues = eig(scaled);
catch
    bound = 0.9/magnitude;
    return;
end
if any(~isfinite(eigenvalues))
    bound = 0.9/magnitude;
    return;
end
imaginary_tolerance = 100*eps*max(1, magnitude);
negative_real = real(eigenvalues) < 0 & ...
                abs(imag(eigenvalues)) <= imaginary_tolerance;
if any(negative_real)
    first_pole = min(-1./real(eigenvalues(negative_real)));
    bound = 0.9*first_pole;
else
    bound = Inf;
end
end
