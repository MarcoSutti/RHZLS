function [problem, x0, reference] = make_problem(name, config)
%MAKE_PROBLEM Build a spectral example for the same RHZLS solver.
%   sphere:    n (60), optional symmetric A and x0.
%   stiefel:   n (60), p (3), optional symmetric A, weights and x0.
    if nargin < 2 || isempty(config), config = struct(); end
    if ~isstruct(config) || ~isscalar(config)
        error('rhzls:InvalidProblem', 'Problem settings must be a scalar structure.');
    end
    name = lower(char(name));
    switch name
        case {'sphere', 'stiefel'}
            n = option(config, 'n', 60);
            if isfield(config, 'A')
                A = config.A;
                n = size(A, 1);
            else
                validate_dimension(n, 'n');
                A = randn(n) / sqrt(n);
                A = (A + A') / 2;
            end
            validate_symmetric(A);
            eigenvalues = sort(eig(A), 'ascend');
            if strcmp(name, 'sphere')
                [M, curve, maxstep] = rhzls.manifold(name, n);
                problem.cost = @(z) z' * (A * z);
                problem.egrad = @(z) 2 * A * z;
                reference = struct('optimum', eigenvalues(1), 'A', A, ...
                    'eigenvalues', eigenvalues, 'description', 'Minimum Rayleigh quotient');
                problem.metrics = @(z) sphere_metrics(z, A, reference.optimum);
            else
                p = option(config, 'p', min(3, n));
                validate_dimension(p, 'p');
                if p > n, error('rhzls:InvalidProblem', 'Stiefel p must not exceed n.'); end
                weights = option(config, 'weights', p:-1:1);
                validateattributes(weights, {'numeric'}, {'real', 'finite', 'vector', 'numel', p});
                weights = weights(:);
                if any(weights <= 0) || any(diff(weights) > 0)
                    error('rhzls:InvalidProblem', 'Brockett weights must be positive and nonincreasing.');
                end
                N = diag(weights);
                [M, curve, maxstep] = rhzls.manifold(name, [n, p]);
                problem.cost = @(z) sum(sum(z .* ((A * z) * N)));
                problem.egrad = @(z) 2 * (A * z) * N;
                optimum = sum(eigenvalues(1:p) .* weights);
                reference = struct('optimum', optimum, 'A', A, 'weights', weights, ...
                    'eigenvalues', eigenvalues, 'description', 'Minimum Brockett function');
                problem.metrics = @(z) stiefel_metrics(z, A, N, optimum);
            end
            x0 = option(config, 'x0', M.rand());

        otherwise
            error('rhzls:UnknownExample', 'Choose stiefel or sphere (or all in run_rhzls).');
    end
    problem.M = M;
    problem.curve = curve;
    problem.maxstep = maxstep;
    problem.name = name;
    problem.reference = reference;
end

function value = option(config, name, fallback)
    if isfield(config, name), value = config.(name); else, value = fallback; end
end

function validate_dimension(value, name)
    validateattributes(value, {'numeric'}, {'real', 'finite', 'scalar', 'integer', 'positive'}, ...
        'rhzls.make_problem', name);
end

function validate_symmetric(A)
    validateattributes(A, {'numeric'}, {'real', 'finite', '2d', 'square'});
    if isempty(A) || norm(A - A', 'fro') > 1e-12 * max(1, norm(A, 'fro'))
        error('rhzls:InvalidProblem', 'A must be a real symmetric, nonempty matrix.');
    end
end

function metrics = sphere_metrics(x, A, optimum)
    metrics = struct('objective_gap', x' * (A * x) - optimum, ...
        'constraint_error', abs(x' * x - 1));
end

function metrics = stiefel_metrics(x, A, N, optimum)
    metrics = struct('objective_gap', sum(sum(x .* ((A * x) * N))) - optimum, ...
        'constraint_error', norm(x' * x - eye(size(x, 2)), 'fro'));
end

