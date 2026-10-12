function [alpha, y, stats] = linesearch(problem, x, d, user, f0, df0)
%LINESEARCH Apply HZ to phi(t)=cost(curve(x,d,t)).
% problem.curve returns [point, exact ambient curve velocity]. It must be
% the SAME function handle as problem.M.retr. A vector transport alone is
% not a valid replacement for that velocity.
% Optional problem.curve_slope(point,egrad,velocity) must compute the same
% pullback derivative, e.g. using a stable Riemannian metric pairing.
if nargin < 4, user = struct(); end
o = rhzls.options(user);
if ~isfield(problem, 'curve') || ~isequal(problem.M.retr, problem.curve)
    error('rhzls:CurveMismatch', ...
        'Set problem.M.retr and problem.curve to the same two-output curve handle.');
end
if isfield(problem, 'maxstep')
    cap = problem.maxstep(x, d);
    if ~isreal(cap) || ~isscalar(cap) || isnan(cap) || cap <= 0
        error('rhzls:MaxStep', 'problem.maxstep must return a positive scalar or Inf.');
    end
    o.alphamax = min(o.alphamax, cap);
end
cached_t = NaN; cached_y = x;
cost_calls = 0; gradient_calls = 0;
if nargin < 6
    [alpha, ~, ~, stats] = rhzls.hz(@evaluate, o);
    origin_evaluations = 1;
else
    [alpha, ~, ~, stats] = rhzls.hz(@evaluate, o, f0, df0);
    origin_evaluations = 0;
end
stats.origin_evaluations = origin_evaluations;
% Oracle attempts can fail in the retraction before either callback is used.
% Keep their budget separate from actual objective and gradient evaluations.
stats.costevals = cost_calls;
stats.gradevals = gradient_calls;
y = x;
if stats.accepted
    if cached_t ~= alpha
        error('rhzls:Internal', 'Accepted point was not the last oracle evaluation.');
    end
    y = cached_y;
end

    function [ft, dft] = evaluate(t)
        try
            [z, velocity] = problem.curve(x, d, t);
        catch err
            if strcmp(err.identifier, 'rhzls:InvalidRetraction') ...
                    || strcmp(err.identifier, 'rhzls:SingularRetraction')
                ft = Inf; dft = NaN; return;
            end
            rethrow(err);
        end
        cost_calls = cost_calls + 1;
        ft = problem.cost(z);
        gradient_calls = gradient_calls + 1;
        egrad = problem.egrad(z);
        if isfield(problem, 'curve_slope')
            dft = problem.curve_slope(z, egrad, velocity);
        else
            dft = rhzls.ambient_inner(egrad, velocity);
        end
        cached_t = t; cached_y = z;
    end
end
