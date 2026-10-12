function report = run_rhzls_tests()
%RUN_RHZLS_TESTS Deterministic regressions; a failed check raises an error.
% Finite differences independently verify the exact retraction velocities.
previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(19, 'twister');
tests = {
    'HZ expands a small first step', @test_expansion;
    'HZ brackets an overshoot', @test_overshoot;
    'Approximate Wolfe at a large objective baseline', @test_approximate;
    'Approximate-only rejects a rounded-cost Wolfe overshoot', @test_approximate_only;
    'HZ U3 contracts a high value with negative slope', @test_u3;
    'Unbounded decrease stops at alphamax', @test_alphamax;
    'Non-descent directions are rejected', @test_not_descent;
    'Oracle budget includes an uncached origin', @test_budget;
    'Nonfinite domain trials are contracted', @test_invalid_domain;
    'Sphere curve velocity and tangent constraint', @() test_geometry('sphere', 9);
    'Stiefel polar velocity and tangent constraint', @() test_geometry('stiefel', [11 3]);
    'Fixed-rank AO velocity and tangent constraint', @() test_geometry('fixedrank', [12 10 3]);
    'Factored ambient pairings agree with dense pairings', @test_pairings;
    'Sphere pullback derivative and accepted step', @() test_pullback('sphere', struct('n', 9));
    'Brockett pullback derivative and accepted step', @() test_pullback('stiefel', struct('n', 11, 'p', 3));
    'Fixed-rank step stays in the smooth retraction domain', @test_rank_domain;
    'Fixed-rank chart exhaustion rejects a non-Wolfe step', @test_rank_domain_exhaustion;
    'Nonnormal fixed-rank core permits steps beyond a norm bound', @test_nonnormal_rank_domain;
    'Mismatched retraction and derivative are rejected', @test_curve_mismatch;
    'Sphere optimization reaches its spectral reference', @test_sphere_convergence;
    'Brockett optimization reaches its spectral reference', @test_stiefel_convergence;
    'Manopt adapter returns a valid step and cache key', @test_manopt_adapter;
    'Sphere comparison shares the problem and initial point', @() test_comparison('sphere');
    'Brockett comparison shares the problem and initial point', @() test_comparison('stiefel');
    'Comparison distances respect optimizer sign and basis ambiguity', @test_comparison_gauges;
    'Default Armijo accepts its first valid trial without halving', @test_default_armijo;
    'RHZLS callback counts distinguish cached and uncached origins', @test_rhzls_evaluation_counts;
    'Invalid retractions consume oracle budget without callback evaluations', @test_invalid_retraction_counts;
    'Armijo callback counts and comparison nf/ng histories', @test_armijo_evaluation_counts;
    'Missing gradient counts require the official default Armijo callback', @test_armijo_count_contract
    };
report = struct('passed', 0, 'failed', 0, 'names', {tests(:, 1)}, 'failures', {{}});
fprintf('RHZLS regression suite: %d tests\n', size(tests, 1));
for k = 1:size(tests, 1)
    try
        tests{k, 2}();
        report.passed = report.passed + 1;
        fprintf('  PASS  %s\n', tests{k, 1});
    catch err
        report.failed = report.failed + 1;
        report.failures{end+1} = struct('name', tests{k, 1}, ...
            'identifier', err.identifier, 'message', err.message); %#ok<AGROW>
        fprintf(2, '  FAIL  %s\n        %s [%s]\n', tests{k, 1}, err.message, err.identifier);
    end
end
fprintf('RHZLS tests: %d passed, %d failed.\n', report.passed, report.failed);
if report.failed > 0
    error('rhzls:TestsFailed', '%d of %d RHZLS regression tests failed.', report.failed, size(tests, 1));
end
end

function test_expansion()
o = rhzls.options(struct('alpha0', .01));
oracle = @(t) quadratic(t, 3);
[a, f, df, stats] = rhzls.hz(oracle, o, 9, -6);
verify(stats.accepted && a > o.alpha0 && stats.evaluations > 1, ...
    'A small initial step should be expanded before acceptance.');
verify_wolfe(a, f, df, 9, -6, o, stats);
end

function test_overshoot()
o = rhzls.options(struct('alpha0', 10));
[a, f, df, stats] = rhzls.hz(@(t) quadratic(t, 1), o, 1, -2);
verify(stats.accepted && a > 0 && a < o.alpha0, ...
    'An overshoot should be bracketed and contracted.');
verify_wolfe(a, f, df, 1, -2, o, stats);
end

function test_approximate()
o = rhzls.options(struct('alpha0', 1));
[f0, df0] = plateau(0);
[a, f, df, stats] = rhzls.hz(@plateau, o, f0, df0);
verify(stats.accepted && stats.approxwolfe && ~stats.wolfe, ...
    'A permitted small objective increase should use approximate Wolfe only.');
verify(f > f0 && f <= stats.phi_limit, 'The approximate-Wolfe objective allowance was not applied.');
verify_wolfe(a, f, df, f0, df0, o, stats);
end

function test_u3()
times = [];
o = rhzls.options(struct('alpha0', 5, 'theta', .5));
[a, f, df, stats] = rhzls.hz(@oracle, o, 0, -1);
verify(numel(times) >= 2 && times(1) == 5 && times(2) == 2.5, ...
    'The high-cost, negative-slope U3 case must bisect before interpolation.');
verify_wolfe(a, f, df, 0, -1, o, stats);
    function [f, df] = oracle(t)
        times(end+1) = t;
        f = -t + 3*t^2 - .4*t^3;
        df = -1 + 6*t - 1.2*t^2;
    end
end

function test_alphamax()
o = rhzls.options(struct('alpha0', .01, 'alphamax', .2));
[a, f, df, stats] = rhzls.hz(@linear_decrease, o, 0, -1);
verify(~stats.accepted && a == 0 && f == 0 && df == -1, ...
    'A search without Wolfe curvature must return its origin.');
verify(strcmp(stats.status, 'alphamax'), 'The finite expansion limit should be reported.');
end

function test_not_descent()
for slope = [0 1]
    [a, f, df, stats] = rhzls.hz(@linear_decrease, struct(), 7, slope);
    verify(~stats.accepted && a == 0 && f == 7 && df == slope ...
        && stats.evaluations == 0 && strcmp(stats.status, 'not_descent'), ...
        'A cached zero or positive initial slope must be rejected without trial evaluations.');
end
end

function test_budget()
o = rhzls.options(struct('alpha0', 10, 'maxevals', 1));
[a, ~, ~, stats] = rhzls.hz(@(t) quadratic(t, 1), o, 1, -2);
verify(~stats.accepted && a == 0 && stats.evaluations == 1 ...
    && strcmp(stats.status, 'max_evaluations'), 'A cached origin permits exactly one trial.');
[a, ~, ~, stats] = rhzls.hz(@(t) quadratic(t, 1), o);
verify(~stats.accepted && a == 0 && stats.evaluations == 1 ...
    && strcmp(stats.status, 'max_evaluations'), 'An uncached origin consumes one evaluation.');
end

function test_invalid_domain()
o = rhzls.options(struct('alpha0', 2));
[a, f, df, stats] = rhzls.hz(@oracle, o, .04, -.4);
verify(stats.accepted && a > 0 && a < .5 && stats.evaluations >= 4, ...
    'Nonfinite domain trials should be contracted to a finite accepted step.');
verify_wolfe(a, f, df, .04, -.4, o, stats);
    function [f, df] = oracle(t)
        if t >= .5
            f = Inf; df = NaN;
        else
            [f, df] = quadratic(t, .2);
        end
    end
end

function test_geometry(name, dimensions)
[M, curve, cap] = rhzls.manifold(name, dimensions);
x = M.rand();
d = M.randvec(x);
d = M.lincomb(x, 1/M.norm(x, d), d);
limit = cap(x, d);
t = min(.2, .3*limit);
times = [0 t];
for t = times
    [y, velocity] = curve(x, d, t);
    Y = rhzls.full_matrix(y);
    V = rhzls.full_matrix(velocity);
    h = min(2e-6, .05*limit);
    if t == 0
        Y1 = rhzls.full_matrix(curve(x, d, h));
        Y2 = rhzls.full_matrix(curve(x, d, 2*h));
        fd = (-3*Y + 4*Y1 - Y2)/(2*h);
        if isstruct(d)
            D = x.U*d.M*x.V' + d.Up*x.V' + x.U*d.Vp';
        else
            D = d;
        end
        relative_close(V, D, 2e-11, 'The derivative at zero must equal the tangent direction.');
    else
        h = min(h, t/4);
        Yp = rhzls.full_matrix(curve(x, d, t+h));
        Ym = rhzls.full_matrix(curve(x, d, t-h));
        fd = (Yp-Ym)/(2*h);
    end
    relative_close(V, fd, 3e-6, 'The analytic curve velocity disagrees with finite differences.');
    relative_close(Y, rhzls.full_matrix(M.retr(x, d, t)), 1e-12, ...
        'The optimized retraction and differentiated curve must agree.');
    verify_constraint(name, y, dimensions);
    switch name
        case 'sphere'
            residual = abs(Y(:)'*V(:));
        case 'stiefel'
            residual = norm(Y'*V + V'*Y, 'fro');
        case 'fixedrank'
            residual = norm(V - y.U*(y.U'*V) - (V*y.V)*y.V' ...
                + y.U*(y.U'*V*y.V)*y.V', 'fro');
    end
    verify(residual <= 1e-10*max(1, norm(V, 'fro')), ...
        'A curve velocity must lie in the tangent space at its returned point.');
end
end

function test_pairings()
G = randn(7, 3); H = randn(5, 3);
v = struct('G', G, 'H', H);
E = sprandn(7, 5, .35);
reference = sum(sum(full(E).*(G*H')));
relative_close(rhzls.ambient_inner(E, v), reference, 1e-12, 'Sparse/factored pairing mismatch.');
eg = struct('U', randn(7, 2), 'S', randn(2, 2), 'V', randn(5, 2));
Efull = rhzls.full_matrix(eg);
reference = sum(sum(Efull.*(G*H')));
relative_close(rhzls.ambient_inner(eg, v), reference, 1e-12, 'Triplet/factored pairing mismatch.');
relative_close(rhzls.ambient_inner(eg, G*H'), reference, 1e-12, 'Triplet/dense pairing mismatch.');
end

function test_pullback(name, config)
[p, x] = rhzls.make_problem(name, config);
gradient = p.M.egrad2rgrad(x, p.egrad(x));
d = p.M.lincomb(x, -1/p.M.norm(x, gradient), gradient);
t = min(.15, .25*p.maxstep(x, d));
[y, v] = p.curve(x, d, t);
analytic = rhzls.ambient_inner(p.egrad(y), v);
h = min(2e-6, t/4);
yp = p.curve(x, d, t+h); ym = p.curve(x, d, t-h);
fd = (p.cost(yp)-p.cost(ym))/(2*h);
relative_close(analytic, fd, 3e-6, 'The pullback derivative disagrees with cost finite differences.');
[~, v0] = p.curve(x, d, 0);
df0 = p.M.inner(x, gradient, d);
relative_close(rhzls.ambient_inner(p.egrad(x), v0), df0, 1e-10, ...
    'The initial pullback slope differs from the Riemannian inner product.');
o = rhzls.options();
f0 = p.cost(x);
[a, z, stats] = rhzls.linesearch(p, x, d, o, f0, df0);
verify_wolfe(a, stats.f, stats.df, f0, df0, o, stats);
relative_close(p.cost(z), stats.f, 1e-12, 'Accepted point and reported cost disagree.');
end

function test_rank_domain()
[M, curve, cap] = rhzls.manifold('fixedrank', [8 7 2]);
x = M.rand();
d = struct('Up', zeros(size(x.U)), 'M', -x.S, 'Vp', zeros(size(x.V)));
X = rhzls.full_matrix(x);
normsq = sum(X(:).^2);
p = struct('M', M, 'curve', curve, 'maxstep', cap, ...
    'cost', @(z) .5*sum(sum(rhzls.full_matrix(z).^2)), ...
    'egrad', @rhzls.full_matrix);
o = rhzls.options(struct('alpha0', 2, 'alphamax', 3));
[a, y, stats] = rhzls.linesearch(p, x, d, o, .5*normsq, -normsq);
verify(stats.accepted && a <= cap(x, d) && a < 1, ...
    'The AO singularity at t=1 must remain outside the search domain.');
relative_close(rhzls.full_matrix(y), (1-a)*X, 1e-12, 'The capped AO trial is incorrect.');
verify_wolfe(a, stats.f, stats.df, .5*normsq, -normsq, o, stats);
end

function test_curve_mismatch()
[p, x] = rhzls.make_problem('sphere', struct('n', 6));
d = p.M.randvec(x);
p.M.retr = @(x, d, t) x; %#ok<INUSD>
expect_error(@() rhzls.linesearch(p, x, d), 'rhzls:CurveMismatch');
end

function test_rank_domain_exhaustion()
saved_rng = rng;
restore_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
[M, curve, cap] = rhzls.manifold('fixedrank', [8 7 2]);
x = M.rand();
d = struct('Up', zeros(size(x.U)), 'M', -x.S, 'Vp', zeros(size(x.V)));
X = rhzls.full_matrix(x);
normsq = sum(X(:).^2);
p = struct('M', M, 'curve', curve, 'maxstep', cap, ...
    'cost', @(z) .5*sum(sum((rhzls.full_matrix(z)+X).^2)), ...
    'egrad', @(z) rhzls.full_matrix(z)+X);
% Along this curve phi(t)=.5*(2-t)^2*||X||^2. With sigma=.1,
% curvature requires t>=1.8, but the first AO pole is t=1.
o = rhzls.options(struct('alpha0', 2, 'alphamax', 3, 'sigma', .1));
[a, y, stats] = rhzls.linesearch(p, x, d, o, 2*normsq, -2*normsq);
verify(cap(x, d) < 1 && ~stats.accepted && a == 0 ...
    && strcmp(stats.status, 'alphamax'), ...
    'Exhausting the smooth rank chart must report failure, without accepting an Armijo-only step.');
relative_close(rhzls.full_matrix(y), X, 1e-12, 'A failed line search must return its input point.');
relative_close(stats.f, 2*normsq, 1e-12, 'A failed line search must report its input cost.');
end

function test_nonnormal_rank_domain()
saved_rng = rng;
restore_rng = onCleanup(@() rng(saved_rng)); %#ok<NASGU>
[M, curve, cap] = rhzls.manifold('fixedrank', [8 7 2]);
x = M.rand(); x.S = eye(2);
d = struct('Up', zeros(size(x.U)), 'M', [0 20; 0 0], 'Vp', zeros(size(x.V)));
X = rhzls.full_matrix(x);
D = x.U*d.M*x.V';
target = X+D;
normsq = sum(D(:).^2);
% S\M is nilpotent: det(S+t*M)=1 for every real t. Its large
% operator norm cannot be used as the true endpoint of the rank chart.
verify(isinf(cap(x, d)), 'A nilpotent AO core has no positive singularity.');
[y, v] = curve(x, d, .5);
relative_close(rhzls.full_matrix(y), X+.5*D, 1e-12, 'A nonnormal AO trial is incorrect.');
relative_close(rhzls.full_matrix(v), D, 1e-12, 'A nonnormal AO velocity is incorrect.');
p = struct('M', M, 'curve', curve, 'maxstep', cap, ...
    'cost', @(z) .5*sum(sum((rhzls.full_matrix(z)-target).^2)), ...
    'egrad', @(z) rhzls.full_matrix(z)-target);
o = rhzls.options(struct('alpha0', 1));
[a, z, stats] = rhzls.linesearch(p, x, d, o, .5*normsq, -normsq);
verify_wolfe(a, stats.f, stats.df, .5*normsq, -normsq, o, stats);
verify(a > .9/norm(x.S\d.M, 2), 'A sufficient-only norm estimate blocked a valid Wolfe step.');
relative_close(rhzls.full_matrix(z), target, 1e-12, 'The nonnormal AO line search missed its quadratic minimum.');
end

function test_sphere_convergence()
A = diag([1 2 4 6 9 12 15 20]);
[p, x0, reference] = rhzls.make_problem('sphere', struct('n', 8, 'A', A));
[x, f, info, details] = rhzls.solve(p, x0, solver_options(1000));
verify_constraint('sphere', x, 8);
verify(abs(f-reference.optimum) < 2e-7, 'Rayleigh minimization did not reach the smallest eigenvalue.');
verify(info(end).gradnorm < 2e-5, 'Rayleigh minimization has not reached stationarity.');
verify_history(info, details);
end

function test_stiefel_convergence()
A = diag([1 2 3 5 7 10 13 17 21 25]);
[p, x0, reference] = rhzls.make_problem('stiefel', struct('n', 10, 'p', 3, 'A', A, 'weights', [3 2 1]));
[x, f, info, details] = rhzls.solve(p, x0, solver_options(2000));
verify_constraint('stiefel', x, [10 3]);
verify(abs(f-reference.optimum) < 2e-6, 'Brockett minimization did not reach its ordered spectral minimum.');
verify(info(end).gradnorm < 1e-4, 'Brockett minimization has not reached stationarity.');
verify_history(info, details);
end

function test_manopt_adapter()
[p, x] = rhzls.make_problem('sphere', struct('n', 8));
g = p.M.egrad2rgrad(x, p.egrad(x));
d = p.M.lincomb(x, -1, g);
df0 = -p.M.norm(x, g)^2;
storedb = StoreDB();
key = storedb.getNewKey();
options = struct('rhzls', struct());
[step, y, newkey, stats] = rhzls.manopt_linesearch(p, x, d, p.cost(x), df0, options, storedb, key);
verify(stats.accepted && step > 0 && ~isequal(newkey, key), 'The Manopt adapter must allocate a key for an accepted point.');
verify_wolfe(step/p.M.norm(x, d), stats.f, stats.df, p.cost(x), df0, rhzls.options(), stats);
relative_close(p.cost(y), stats.f, 1e-12, 'Manopt adapter returned an inconsistent point.');
solver = struct('linesearch', @rhzls.manopt_linesearch, 'rhzls', struct(), ...
    'maxiter', 600, 'tolgradnorm', 1e-6, 'minstepsize', 0, 'verbosity', 0);
[z, f, info] = steepestdescent(p, x, solver);
verify_constraint('sphere', z, 8);
verify(abs(f-p.reference.optimum) < 2e-7 && info(end).gradnorm < 2e-5, ...
    'Manopt steepestdescent did not converge with the RHZLS adapter.');
for k = 2:numel(info)
    s = info(k).linesearch;
    verify(s.accepted && (s.wolfe || s.approxwolfe), 'Manopt accepted an invalid RHZLS iteration.');
end
end

function test_comparison(name)
switch name
    case 'sphere'
        config = struct('n', 6, 'A', diag([1 2 3 5 8 11]));
        dimensions = 6;
    case 'stiefel'
        config = struct('n', 7, 'p', 2, 'A', diag([1 2 3 4 6 9 12]), 'weights', [2 1]);
        dimensions = [7 2];
end
[p, x0] = rhzls.make_problem(name, config);
o = struct('maxiter', 12, 'tolgradnorm', 0, 'minstepsize', 0, 'verbosity', 0);
result = rhzls.compare(p, x0, o);
initial_cost = p.cost(x0);
g0 = p.M.egrad2rgrad(x0, p.egrad(x0));
initial_gradnorm = p.M.norm(x0, g0);
verify(initial_gradnorm > 0, 'The comparison fixture must have a nonstationary initial point.');
relative_close(rhzls.full_matrix(result.x0), rhzls.full_matrix(x0), 1e-14, ...
    'The common comparison point differs from the supplied point.');
relative_close(result.reference.gradient_scale, initial_gradnorm, 1e-12, ...
    'The relative-gradient normalization must use the common initial norm.');
methods = {'rhzls', 'manopt'};
callbacks = {'rhzls.manopt_linesearch', 'linesearch'};
for j = 1:numel(methods)
    run = result.(methods{j});
    history = run.info;
    verify(numel(history) > 1 && all([history.iter] == 0:numel(history)-1) ...
        && history(end).iter <= o.maxiter, 'The comparison history has invalid iteration indices.');
    verify(run.details.options.maxiter == o.maxiter ...
        && run.details.options.tolgradnorm == o.tolgradnorm ...
        && run.details.options.minstepsize == o.minstepsize, ...
        'Both methods must use the same stopping criteria.');
    verify(strcmp(func2str(run.details.options.linesearch), callbacks{j}), ...
        'The comparison selected a different line-search callback.');
    relative_close(history(1).cost, initial_cost, 1e-12, ...
        'A method did not start at the common initial cost.');
    relative_close(history(1).gradnorm, initial_gradnorm, 1e-12, ...
        'A method did not start at the common initial gradient.');
    expected_objective = abs([history.cost] - p.reference.optimum);
    relative_close([history.objective_error], expected_objective, 1e-12, ...
        'Absolute objective-error history does not match its recorded costs and optimum.');
    relative_close([history.gradient_error], [history.gradnorm] / initial_gradnorm, 1e-12, ...
        'Gradient-error history uses an inconsistent normalization.');
    verify(all(isfinite([history.solution_error])) && all([history.solution_error] >= 0), ...
        'The comparison solution errors must be finite and nonnegative.');
    expected_solution = independent_solution_error(name, run.x, p.reference);
    relative_close(history(end).solution_error, expected_solution, 1e-11, ...
        'The final solution-error history does not describe the final point.');
    relative_close(history(end).cost, run.cost, 1e-12, ...
        'The final comparison cost disagrees with its history.');
    relative_close(p.cost(run.x), run.cost, 1e-12, ...
        'The final comparison cost disagrees with its returned point.');
    verify_constraint(name, run.x, dimensions);
end
relative_close(result.rhzls.info(1).solution_error, result.manopt.info(1).solution_error, ...
    1e-14, 'The methods use different solution-reference normalizations.');
% Do not select a baseline callback in this independent replay: the solver
% must choose its own official default, and comparison bookkeeping must
% preserve its trajectory and retraction.
[xa, fa, ia, replay_options] = steepestdescent(p, x0, o);
verify(isequal(replay_options.linesearch, @linesearch), ...
    'The independent solver replay did not select official default Armijo.');
verify(numel(ia) == numel(result.manopt.info), 'The default-Armijo replay has a different history length.');
relative_close([ia.cost], [result.manopt.info.cost], 1e-12, ...
    'The comparison baseline differs from an independent stock-Armijo run.');
relative_close([ia.gradnorm], [result.manopt.info.gradnorm], 1e-12, ...
    'The comparison baseline has different gradients from an independent replay.');
relative_close(rhzls.full_matrix(xa), rhzls.full_matrix(result.manopt.x), 1e-12, ...
    'The comparison baseline and its independent replay return different points.');
relative_close(fa, result.manopt.cost, 1e-12, 'The independent replay has a different final cost.');
if strcmp(name, 'sphere')
    [xh, fh, ih] = rhzls.solve(p, x0, o);
    verify(numel(ih) == numel(result.rhzls.info), ...
        'The Manopt RHZLS adapter changed the standalone iteration count.');
    relative_close([ih.cost], [result.rhzls.info.cost], 1e-11, ...
        'The adapter warm start changed the standalone RHZLS cost trajectory.');
    relative_close([ih.gradnorm], [result.rhzls.info.gradnorm], 1e-11, ...
        'The adapter warm start changed the standalone RHZLS gradient trajectory.');
    relative_close(rhzls.full_matrix(xh), rhzls.full_matrix(result.rhzls.x), 1e-11, ...
        'The adapter and standalone RHZLS return different points.');
    relative_close(fh, result.rhzls.cost, 1e-11, ...
        'The adapter and standalone RHZLS return different final costs.');
end
hz_history = result.rhzls.info;
hz_options = rhzls.options(result.rhzls.details.options.rhzls);
for k = 2:numel(hz_history)
    s = hz_history(k).linesearch;
    verify_wolfe(hz_history(k).alpha, s.f, s.df, s.f0, s.df0, hz_options, s);
end
end

function value = independent_solution_error(name, x, reference)
switch name
    case 'sphere'
        [V, D] = eig(reference.A);
        [~, first] = min(diag(D));
        target = V(:, first);
        value = min(norm(x-target), norm(x+target));
    case 'stiefel'
        [V, D] = eig(reference.A);
        [~, order] = sort(diag(D), 'ascend');
        target = V(:, order(1:size(x, 2)));
        squared_error = 0;
        for k = 1:size(x, 2)
            squared_error = squared_error + min(norm(x(:, k)-target(:, k)), ...
                norm(x(:, k)+target(:, k)))^2;
        end
        value = sqrt(squared_error/size(x, 2));
end
end

function test_comparison_gauges()
% A minimizer and a sign-equivalent representation must have zero error.
E = eye(4);
[p, x0, ref] = rhzls.make_problem('sphere', struct('A', diag([1 2 3 4]), 'x0', -E(:, 1)));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
[p, x0, ref] = rhzls.make_problem('sphere', struct('A', diag([1 1 3 4]), ...
    'x0', -.6*E(:, 1)-.8*E(:, 2)));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
other = metric(E(:, 3));
verify(other.solution_error > .1, 'A vector outside the minimum eigenspace cannot have zero error.');
[p, x0, ref] = rhzls.make_problem('stiefel', struct('A', diag([1 2 3 4]), ...
    'p', 2, 'weights', [2 1], 'x0', [-E(:, 1), E(:, 2)]));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
other = metric(E(:, [2 1]));
verify(other.solution_error > .1, 'Unequal Brockett weights do not permit arbitrary column permutations.');
rotation = [cos(.6), -sin(.6); sin(.6), cos(.6)];
rotated = E(:, 1:2)*rotation;
[p, x0, ref] = rhzls.make_problem('stiefel', struct('A', diag([1 2 3 4]), ...
    'p', 2, 'weights', [1 1], 'x0', rotated));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
[p, x0, ref] = rhzls.make_problem('stiefel', struct('A', diag([1 1 3 4]), ...
    'p', 2, 'weights', [2 1], 'x0', rotated));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
[p, x0, ref] = rhzls.make_problem('stiefel', struct('A', diag([1 2 2 4]), ...
    'p', 2, 'weights', [2 1], 'x0', E(:, [1 3])));
metric = rhzls.comparison_metrics(p, ref, x0);
verify_zero_metric(metric(x0));
end

function verify_zero_metric(metric)
verify(metric.objective_error < 2e-12 && metric.solution_error < 2e-12, ...
    'An optimizer-equivalent representation has a nonzero comparison error.');
end

function test_approximate_only()
oracle = @(a) deal(1e20 + .5*(a-1)^2, a-1);
[a,~,~,stats] = rhzls.hz(oracle,struct('alpha0',10,'acceptance','approximate'));
verify(stats.accepted && stats.approxwolfe, 'Approximate-only must satisfy both derivative bounds.');
verify(abs(a-1) < 1e-12, 'A rounded objective must not accept the overshoot at alpha=10.');
end

function test_default_armijo()
[p, x] = rhzls.make_problem('sphere', struct('A', diag([1 2]), 'x0', [1; 1]/sqrt(2)));
g = p.M.egrad2rgrad(x, p.egrad(x));
d = p.M.lincomb(x, -1, g);
normsq = p.M.norm(x, g)^2;
storedb = StoreDB(); key = storedb.getNewKey();
[step, y, newkey, stats] = linesearch(p, x, d, p.cost(x), ...
    -normsq, struct(), storedb, key);
verify(stats.costevals == 1 && ~isequal(newkey, key), ...
    'A valid first Armijo trial must be accepted without an extra trial.');
relative_close(stats.alpha, 1, 1e-12, 'Official Armijo must not halve an already valid first trial.');
relative_close(step, 1, 1e-12, 'The default first displacement has the wrong norm.');
relative_close(y, [1; 0], 1e-12, 'The default first trial did not retract to the expected minimizer.');
relative_close(p.cost(y), 1, 1e-12, 'The default first trial has the wrong objective value.');
end

function test_rhzls_evaluation_counts()
% Observe actual callback invocations independently of the oracle budget.
cost_points = []; gradient_points = [];
curvefun = @curve;
M = struct('retr', curvefun, 'norm', @(z, v) abs(v));
p = struct('M', M, 'curve', curvefun, 'cost', @cost, 'egrad', @egrad);
o = struct('alpha0', .01);
for cached = [true false]
    cost_points = []; gradient_points = [];
    if cached
        [alpha, ~, stats] = rhzls.linesearch(p, 0, 1, o, 1, -2);
    else
        [alpha, ~, stats] = rhzls.linesearch(p, 0, 1, o);
    end
    verify(stats.accepted && alpha > 0 && numel(cost_points) > 1, ...
        'The accounting fixture must execute multiple genuine trial callbacks.');
    verify(stats.costevals == numel(cost_points) ...
        && stats.gradevals == numel(gradient_points), ...
        'RHZLS evaluation counts disagree with independently observed callback calls.');
    verify(sum(cost_points == 0) == ~cached && sum(gradient_points == 0) == ~cached ...
        && stats.origin_evaluations == ~cached, ...
        'An uncached origin must be counted once; a supplied origin must not call the callbacks.');
end
% The adapter must preserve the counts rather than replace them by a budget.
cost_points = []; gradient_points = [];
storedb = StoreDB(); key = storedb.getNewKey();
[step, ~, newkey, stats] = rhzls.manopt_linesearch(p, 0, 1, 1, -2, ...
    struct('rhzls', o), storedb, key);
verify(step > 0 && ~isequal(key, newkey) ...
    && stats.costevals == numel(cost_points) ...
    && stats.gradevals == numel(gradient_points) ...
    && all(cost_points > 0) && all(gradient_points > 0), ...
    'The Manopt adapter lost actual evaluation counts or reevaluated its cached origin.');

    function [z, velocity] = curve(x, d, t)
        z = x + t*d; velocity = d;
    end
    function value = cost(z)
        cost_points(end+1) = z;
        value = (z-1)^2;
    end
    function value = egrad(z)
        gradient_points(end+1) = z;
        value = 2*(z-1);
    end
end

function test_invalid_retraction_counts()
% A rejected curve evaluation is not an evaluation of the user's objective.
cost_calls = 0; gradient_calls = 0; invalid_trials = 0;
curvefun = @curve;
M = struct('retr', curvefun);
p = struct('M', M, 'curve', curvefun, 'cost', @cost, 'egrad', @egrad);
[alpha, ~, stats] = rhzls.linesearch(p, 0, 1, struct('alpha0', 2), .04, -.4);
verify(stats.accepted && alpha > 0 && alpha < .5 && invalid_trials > 0, ...
    'The invalid-retraction fixture must recover from a rejected trial.');
verify(stats.costevals == cost_calls && stats.gradevals == gradient_calls ...
    && cost_calls > 0 && gradient_calls > 0, ...
    'Actual counts must match callback calls after invalid trials.');
verify(stats.evaluations > stats.costevals ...
    && stats.evaluations == cost_calls + invalid_trials, ...
    'Invalid retractions must consume oracle budget without fictitious cost evaluations.');

    function [z, velocity] = curve(x, d, t)
        if t >= .5
            invalid_trials = invalid_trials + 1;
            error('rhzls:InvalidRetraction', 'Synthetic curve domain ends at t=.5.');
        end
        z = x + t*d; velocity = d;
    end
    function value = cost(z)
        verify(z < .5, 'The cost callback must never receive an invalid trial.');
        cost_calls = cost_calls + 1;
        value = (z-.2)^2;
    end
    function value = egrad(z)
        verify(z < .5, 'The gradient callback must never receive an invalid trial.');
        gradient_calls = gradient_calls + 1;
        value = 2*(z-.2);
    end
end

function test_armijo_evaluation_counts()
[p, x] = rhzls.make_problem('sphere', struct('A', diag([1 2]), 'x0', [1; 1]/sqrt(2)));
base_cost = p.cost; base_egrad = p.egrad;
g = p.M.egrad2rgrad(x, base_egrad(x));
d = p.M.lincomb(x, -1, g);
f0 = base_cost(x); df0 = -p.M.norm(x, g)^2;
cost_calls = 0; gradient_calls = 0;
p.cost = @cost; p.egrad = @egrad;
storedb = StoreDB(); key = storedb.getNewKey();
[step, ~, ~, stats] = linesearch(p, x, d, f0, df0, ...
    struct('ls_initial_stepsize', 1000), storedb, key);
verify(step > 0 && cost_calls > 1, 'The Armijo fixture must backtrack through rejected trials.');
verify(stats.costevals == cost_calls && gradient_calls == 0, ...
    'Armijo must count rejected cost trials and report zero line-search gradient evaluations.');
if isfield(stats, 'gradevals')
    verify(stats.gradevals == gradient_calls, ...
        'Any gradient count returned by official Armijo must agree with observed callbacks.');
end

% Verify the public per-iteration fields, including their k=0 convention.
[problem, x0] = rhzls.make_problem('sphere', struct('A', diag([1 2 5])));
result = rhzls.compare(problem, x0, struct('maxiter', 3, 'tolgradnorm', 0, ...
    'minstepsize', 0, 'verbosity', 0));
for method = {'manopt', 'rhzls'}
    run = result.(method{1});
    history = run.info;
    verify(isfield(history, 'nf') && isfield(history, 'ng') ...
        && history(1).nf == 0 && history(1).ng == 0, ...
        'Comparison histories require nf/ng fields with zero line-search work at iteration zero.');
    for k = 2:numel(history)
        if isfield(history(k).linesearch, 'gradevals')
            expected_ng = history(k).linesearch.gradevals;
        else
            verify(isequal(run.details.options.linesearch, @linesearch), ...
                'Only the verified official cost-only default may omit gradient counts.');
            expected_ng = 0;
        end
        verify(history(k).nf == history(k).linesearch.costevals ...
            && history(k).ng == expected_ng ...
            && history(k).nf >= 1 && history(k).nf == fix(history(k).nf) ...
            && history(k).ng >= 0 && history(k).ng == fix(history(k).ng), ...
            'Comparison nf/ng must preserve the actual integer counts for that line search.');
    end
end
verify(all([result.manopt.info.ng] == 0), ...
    'Outer-solver gradients must not appear in the Armijo line-search ng history.');

    function value = cost(z)
        cost_calls = cost_calls + 1;
        value = base_cost(z);
    end
    function value = egrad(z)
        gradient_calls = gradient_calls + 1;
        value = base_egrad(z);
    end
end

function test_armijo_count_contract()
% A missing gradient field cannot be interpreted as zero for arbitrary
% callbacks. The official Armijo search is independently observed above.
raw = struct('iter', {0, 1}, 'linesearch', {[], struct('costevals', 3)});
expect_error(@() rhzls.add_evaluation_counts(raw), 'rhzls:MissingEvaluationCounts');
expect_error(@() rhzls.add_evaluation_counts(raw, @rhzls.manopt_linesearch), ...
    'rhzls:MissingEvaluationCounts');
counted = rhzls.add_evaluation_counts(raw, @linesearch);
verify(isequal([counted.nf], [0 3]) && isequal([counted.ng], [0 0]), ...
    'The official Armijo contract must preserve trial counts and exclude outer gradients.');
verify(~isfield(counted(2).linesearch, 'gradevals'), ...
    'Evaluation bookkeeping must not rewrite the upstream line-search statistics.');
raw(2).linesearch.gradevals = 2;
counted = rhzls.add_evaluation_counts(raw, @rhzls.manopt_linesearch);
verify(counted(2).ng == 2, 'Explicit callback gradient counts must be preserved.');
end

function o = solver_options(maxiter)
o = struct('maxiter', maxiter, 'tolgradnorm', 1e-7, 'verbosity', 0);
end

function verify_history(info, details)
verify(numel(info) > 1 && ~startsWith(details.status, 'linesearch_'), ...
    'Optimization stopped because its line search failed.');
for k = 2:numel(info)
    s = info(k).linesearch;
    o = rhzls.options(details.options.linesearch);
    verify_wolfe(info(k).alpha, s.f, s.df, s.f0, s.df0, o, s);
end
end

function verify_constraint(name, y, dimensions)
Y = rhzls.full_matrix(y);
switch name
    case 'sphere'
        verify(abs(norm(Y(:))-1) <= 2e-11, 'The point is not on the unit sphere.');
    case 'stiefel'
        verify(norm(Y'*Y-eye(dimensions(2)), 'fro') <= 2e-10, 'The point violates the Stiefel constraint.');
    case 'fixedrank'
        r = dimensions(3);
        verify(norm(y.U'*y.U-eye(r), 'fro') <= 2e-10 ...
            && norm(y.V'*y.V-eye(r), 'fro') <= 2e-10, 'Fixed-rank singular factors are not orthonormal.');
        s = svd(Y);
        verify(s(r) > 1e-10*max(1, s(1)) ...
            && (numel(s) == r || s(r+1) <= 2e-10*max(1, s(1))), ...
            'The point does not have its prescribed rank.');
end
end

function verify_wolfe(a, f, df, f0, df0, o, stats)
% Compute the mathematical acceptance tests independently of stats flags.
limit = f0 + o.epsilon*max(abs(f0), o.function_scale);
ordinary = f <= f0 + o.delta*a*df0 && df >= o.sigma*df0;
approximate = f <= limit && df >= o.sigma*df0 && df <= (2*o.delta-1)*df0;
verify(stats.accepted && a > 0 && isfinite(f) && isfinite(df) ...
    && (ordinary || approximate), 'An accepted step violates both Wolfe and approximate Wolfe.');
verify(stats.wolfe == ordinary && stats.approxwolfe == approximate, 'Acceptance flags disagree with independently computed inequalities.');
relative_close(stats.phi_limit, limit, 2e-15, 'The approximate-Wolfe objective threshold is incorrect.');
end

function [f, df] = quadratic(t, minimizer)
f = (t-minimizer)^2;
df = 2*(t-minimizer);
end

function [f, df] = plateau(t)
f = 1e20 + 1e12*(-t+4*t^2-2.6*t^3);
df = 1e12*(-1+8*t-7.8*t^2);
end

function [f, df] = linear_decrease(t)
f = -t;
df = -1;
end

function relative_close(actual, expected, tolerance, message)
error_norm = norm(actual(:)-expected(:));
scale = max([1 norm(actual(:)) norm(expected(:))]);
verify(error_norm <= tolerance*scale, ...
    sprintf('%s Relative discrepancy: %.3g (limit %.3g).', message, error_norm/scale, tolerance));
end

function expect_error(action, identifier)
try
    action();
catch err
    verify(strcmp(err.identifier, identifier), sprintf('Expected %s but received %s.', identifier, err.identifier));
    return;
end
error('rhzls:AssertionFailed', 'Expected error %s was not raised.', identifier);
end

function verify(condition, message)
if ~condition
    error('rhzls:AssertionFailed', '%s', message);
end
end
