function [x, cost, info, details] = solve(problem, x, user)
%SOLVE Riemannian steepest descent with the common RHZLS line search.
% maxiter, tolgradnorm, minstepsize, verbosity, and linesearch (struct).
if nargin < 3, user = struct(); end
o = rhzls.solver_options(user);
lsoptions = rhzls.options(o.linesearch);
if nargin < 2 || isempty(x), x = problem.M.rand(); end
clock = tic;
cost = problem.cost(x);
grad = problem.M.egrad2rgrad(x, problem.egrad(x));
gn = problem.M.norm(x, grad);
validate_origin(cost, gn);
record = struct('iter', 0, 'cost', cost, 'gradnorm', gn, ...
    'alpha', 0, 'stepsize', 0, 'ls_evaluations', 0, ...
    'ls_status', 'initial', 'time', toc(clock), 'linesearch', struct(), 'nf', 0, 'ng', 0);
info = record;
details = struct('status', 'max_iterations', 'iterations', 0, ...
    'evaluations', 1, 'costevals', 0, 'gradevals', 0, ...
    'failed_linesearch', struct(), 'options', o);
nextalpha = lsoptions.alpha0;
for k = 1:o.maxiter
    if o.tolgradnorm > 0 && gn <= o.tolgradnorm
        details.status = 'converged'; break;
    end
    d = problem.M.lincomb(x, -1, grad);
    current = lsoptions; current.alpha0 = min(nextalpha, current.alphamax);
    [alpha, y, ls] = rhzls.linesearch(problem, x, d, current, cost, -gn^2);
    details.evaluations = details.evaluations + ls.evaluations;
    details.costevals = details.costevals + ls.costevals;
    details.gradevals = details.gradevals + ls.gradevals;
    if ~ls.accepted
        details.status = ['linesearch_' ls.status];
        details.failed_linesearch = ls; break;
    end
    step = alpha*gn;
    x = y; cost = ls.f;
    grad = problem.M.egrad2rgrad(x, problem.egrad(x));
    gn = problem.M.norm(x, grad);
    validate_origin(cost, gn);
    record = struct('iter', k, 'cost', cost, 'gradnorm', gn, ...
        'alpha', alpha, 'stepsize', step, 'ls_evaluations', ls.evaluations, ...
        'ls_status', ls.status, 'time', toc(clock), 'linesearch', ls, ...
        'nf', ls.costevals, 'ng', ls.gradevals);
    info(end+1) = record; %#ok<AGROW>
    details.iterations = k;
    nextalpha = min(2*alpha, lsoptions.alphamax);
    if o.verbosity >= 2
        fprintf('%4d  cost=%+.12e  |grad|=%.3e  alpha=%.3e  %s\n', k, cost, gn, alpha, ls.status);
    end
    if o.tolgradnorm > 0 && gn <= o.tolgradnorm
        details.status = 'converged'; break;
    elseif step <= o.minstepsize
        details.status = 'min_stepsize'; break;
    end
end
if o.tolgradnorm > 0 && gn <= o.tolgradnorm, details.status = 'converged'; end
if o.verbosity >= 1
    fprintf('RHZLS: %s, %d iterations, cost=%+.12e, |grad|=%.3e.\n', ...
        details.status, details.iterations, cost, gn);
end
end

function validate_origin(f, gn)
if ~isreal(f) || ~isscalar(f) || ~isfinite(f) ...
        || ~isreal(gn) || ~isscalar(gn) || ~isfinite(gn)
    error('rhzls:InvalidObjective', 'Cost and gradient norm must be finite real scalars.');
end
end
