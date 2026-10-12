function result = compare(problem, x0, user)
%COMPARE Compare stock Manopt Armijo and RHZLS with the SAME outer solver.
% Both runs use Manopt steepestdescent, the identical problem, x0 and M.retr,
% and common stopping criteria. Only options.linesearch is different.
% Returns result.rhzls and result.manopt, each with x,cost,info,details.
% Each history includes objective_error, gradient_error, and solution_error.
if nargin < 3, user = struct(); end
o = rhzls.solver_options(user);
if nargin < 2 || isempty(x0), x0 = problem.M.rand(); end
if ~isfield(problem, 'reference')
    error('rhzls:MissingReference', 'Comparisons require problem.reference (see rhzls.make_problem).');
end
if ~isfield(problem, 'curve') || ~isequal(problem.M.retr, problem.curve)
    error('rhzls:CurveMismatch', 'The comparison requires matching problem.M.retr and problem.curve.');
end
[metricfun, reference] = rhzls.comparison_metrics(problem, problem.reference, x0);
g0 = problem.M.egrad2rgrad(x0, problem.egrad(x0));
gradient_scale = problem.M.norm(x0, g0);
reference.gradient_scale = gradient_scale;
common = rmfield(o, 'linesearch');
common.statsfun = @record_metrics;

armijo_options = common;
% Leave linesearch unspecified: the official solver selects its own default.
[xa, fa, ia, oa] = steepestdescent(problem, x0, armijo_options);
assert(isequal(oa.linesearch, @linesearch), ...
    'The benchmark baseline must use the official default Armijo search.');
hz_options = common;
hz_options.linesearch = @rhzls.manopt_linesearch;
hz_options.rhzls = o.linesearch;
[xh, fh, ih, oh] = steepestdescent(problem, x0, hz_options);

result = struct('name', problem.name, 'x0', x0, 'reference', reference, 'common_options', o, ...
    'manopt', package_run(xa, fa, ia, oa, 'Manopt default (Armijo)'), ...
    'rhzls', package_run(xh, fh, ih, oh, 'RHZLS'));

    function stats = record_metrics(~, point, stats)
        metrics = metricfun(point, stats.cost);
        stats.objective_error = metrics.objective_error;
        stats.solution_error = metrics.solution_error;
        if gradient_scale > 0
            stats.gradient_error = stats.gradnorm / gradient_scale;
        else
            % The common starting point is stationary; retain absolute norm.
            stats.gradient_error = stats.gradnorm;
        end
    end
end

function run = package_run(x, cost, info, options, method)
evaluations = 1;
for k = 1:numel(info)
    info(k).alpha = 0;
    info(k).ls_evaluations = 0;
    info(k).ls_status = 'initial';
    if k > 1
        stats = info(k).linesearch;
        info(k).alpha = stats.alpha;
        if isfield(stats, 'evaluations')
            info(k).ls_evaluations = stats.evaluations;
            info(k).ls_status = stats.status;
        else
            info(k).ls_evaluations = stats.costevals;
            if stats.alpha > 0
                info(k).ls_status = 'manopt_armijo';
            else
                info(k).ls_status = 'manopt_rejected';
            end
        end
        evaluations = evaluations + info(k).ls_evaluations;
    end
end
info = rhzls.add_evaluation_counts(info, options.linesearch);
last = info(end);
if last.gradnorm < options.tolgradnorm
    status = 'converged';
elseif last.iter >= options.maxiter
    status = 'max_iterations';
elseif last.stepsize < options.minstepsize
    status = 'min_stepsize';
else
    status = 'stopped';
end
details = struct('status', status, 'iterations', last.iter, ...
    'evaluations', evaluations, 'costevals', sum([info.nf]), ...
    'gradevals', sum([info.ng]), 'options', options, 'method', method);
run = struct('x', x, 'cost', cost, 'info', info, 'details', details);
end
