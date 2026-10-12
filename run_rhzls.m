function [x, cost, info, details] = run_rhzls(example, opts)
%RUN_RHZLS Run an educational RHZLS example with a shared solver.
%   run_rhzls                 Brockett function on the Stiefel manifold.
%   run_rhzls('sphere')       Rayleigh quotient on the unit sphere.
%   run_rhzls('all')          Run both spectral examples; outputs have example-name fields.
%
%   opts.seed    Reproducible random seed (default 10).
%   opts.plot    Draw comparison histories, with normalized errors (default true).
%   opts.compare Compare with stock Manopt Armijo (default true).
%   opts.problem Example settings, passed to rhzls.make_problem.
%   opts.solver  Solver settings (sphere maxiter defaults to 250).
%   For 'all', opts.problem may have sphere/stiefel substructures.
%
%   Example:
%     opts = struct('plot', false, 'problem', struct('n', 100));
%     [x, f, info, details] = run_rhzls('sphere', opts);

    if nargin < 1 || isempty(example), example = 'stiefel'; end
    if nargin < 2 || isempty(opts), opts = struct(); end
    if ~isstruct(opts) || ~isscalar(opts)
        error('rhzls:InvalidOptions', 'opts must be a scalar structure.');
    end
    startup_rhzls();
    example = lower(char(example));
    seed = option(opts, 'seed', 10);
    draw_plot = option(opts, 'plot', true);
    do_compare = option(opts, 'compare', true);
    validateattributes(seed, {'numeric'}, {'real', 'finite', 'scalar', 'integer', 'nonnegative'});
    validateattributes(draw_plot, {'logical', 'numeric'}, {'scalar'});
    validateattributes(do_compare, {'logical', 'numeric'}, {'scalar'});
    problem_opts = option(opts, 'problem', struct());
    solver_opts = option(opts, 'solver', struct());
    if ~isstruct(problem_opts) || ~isscalar(problem_opts) || ...
            ~isstruct(solver_opts) || ~isscalar(solver_opts)
        error('rhzls:InvalidOptions', 'opts.problem and opts.solver must be scalar structures.');
    end

    if strcmp(example, 'all')
        names = {'stiefel', 'sphere'};
        x = struct(); cost = struct(); info = struct(); details = struct();
        for k = 1:numel(names)
            local_opts = opts;
            if isfield(problem_opts, names{k})
                local_opts.problem = problem_opts.(names{k});
            end
            [x.(names{k}), cost.(names{k}), info.(names{k}), details.(names{k})] = ...
                run_rhzls(names{k}, local_opts);
        end
        return;
    end

    if strcmp(example, 'sphere') && ~isfield(solver_opts, 'maxiter')
        solver_opts.maxiter = 250;
    end
    % Preserve the caller's random stream after generating reproducible data.
    rng_state = rng;
    restore_rng = onCleanup(@() rng(rng_state)); %#ok<NASGU>
    rng(seed, 'twister');
    [problem, x0, reference] = rhzls.make_problem(example, problem_opts);
    if do_compare
        comparison = rhzls.compare(problem, x0, solver_opts);
        x = comparison.rhzls.x;
        cost = comparison.rhzls.cost;
        info = comparison.rhzls.info;
        details = comparison.rhzls.details;
        details.comparison = comparison;
    else
        [x, cost, info, details] = rhzls.solve(problem, x0, solver_opts);
    end
    details.example = example;
    details.reference = reference;
    details.metrics = problem.metrics(x);
    details.problem = problem;
    details.x0 = x0;
    print_result(example, cost, details);
    if do_compare
        base = comparison.manopt;
        fprintf('  Manopt default (Armijo): %s, %d iterations, |grad|=%.3e\n', ...
            base.details.status, base.details.iterations, base.info(end).gradnorm);
    end
    if draw_plot
        if do_compare
            rhzls.plot_comparison(comparison);
        else
            plot_result(example, info, reference, details.options.tolgradnorm);
        end
    end
end

function value = option(opts, name, fallback)
    if isfield(opts, name), value = opts.(name); else, value = fallback; end
end

function print_result(name, cost, details)
    fprintf('\nRHZLS / %s: cost = %.12g, status = %s\n', name, cost, details.status);
    if isfield(details.metrics, 'objective_gap')
        fprintf('  objective gap = %.3e, constraint error = %.3e\n', ...
            details.metrics.objective_gap, details.metrics.constraint_error);
    end
end

function plot_result(name, info, reference, tolerance)
    if isempty(info), return; end
    rhzls.plot_style(name, numel(info));
    fig = figure('Name', ['RHZLS / ' name], 'NumberTitle', 'off', ...
        'Color', 'w', 'WindowStyle', 'normal', 'Units', 'pixels', 'Position', rhzls.figure_position());
    layout = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    axes_pair = [nexttile(layout, 1), nexttile(layout, 2)];
    rhzls.plot_history(name, info, reference, axes_pair, tolerance);
    title(layout, ['RHZLS / ' name], 'Interpreter', 'latex');
end
