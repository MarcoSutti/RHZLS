function [metricfun, metadata] = comparison_metrics(problem, reference, x0)
%COMPARISON_METRICS Common reference errors for two line-search histories.
%   [metricfun, metadata] = rhzls.comparison_metrics(problem, reference, x0)
%   precomputes the reference solution. metricfun(x)
%   returns objective_error and solution_error. metricfun(x, f) accepts an
%   already computed cost, avoiding an additional objective evaluation.
%
%   The solver supplies gradient_error = gradnorm / initial_gradnorm.
%   objective_error is the absolute gap abs(f - f_star), without scaling.
%   Sphere and Stiefel distances respect the nonuniqueness of minimizers.

    if nargin < 2 || isempty(reference), reference = problem.reference; end
    if nargin < 3, error('rhzls:MissingInitialPoint', 'Supply the common initial point x0.'); end

    f0 = problem.cost(x0);
    known_optimum = isfield(reference, 'optimum') && ...
        isscalar(reference.optimum) && isfinite(reference.optimum);
    metadata = struct('objective_known', known_optimum, ...
        'initial_cost', f0, 'optimum', NaN, 'objective_scale', NaN, ...
        'objective_label', '', 'objective_description', '', ...
        'solution_label', '', 'solution_description', '', ...
        'solution_reference_is_optimum', true);
    if known_optimum
        optimum = reference.optimum;
        metadata.optimum = optimum;
        objective_scale = 1;
        metadata.objective_label = '$|f_k-f_\star|$';
        metadata.objective_description = 'Absolute objective gap, without normalization.';
        metadata.objective_scale = objective_scale;
    else
        optimum = NaN;
        objective_scale = NaN;
        metadata.objective_description = 'Optimal objective unknown; objective_error is NaN.';
    end
    metadata.direct_objective_error = known_optimum && isfield(problem, 'objective_error');

    switch lower(problem.name)
        case 'sphere'
            [vectors, values, eigen_tolerance] = sorted_eigensystem(reference.A);
            minimum_space = vectors(:, values - values(1) <= eigen_tolerance);
            solutionfun = @(x) sphere_distance(x, minimum_space);
            metadata.solution_label = '$d(x_k,\mathcal{X}_\star)$';
            metadata.solution_description = ['Euclidean distance to the nearest unit vector in the minimum ', ...
                'eigenspace; includes the sign ambiguity.'];
            metadata.eigenvalue_tolerance = eigen_tolerance;

        case 'stiefel'
            [vectors, values, eigen_tolerance] = sorted_eigensystem(reference.A);
            weights = reference.weights(:);
            p = numel(weights);
            weight_tolerance = 32 * eps(max(1, max(abs(weights))));
            weight_blocks = equal_value_blocks(weights, weight_tolerance);
            eigen_blocks = equal_value_blocks(values, eigen_tolerance);
            crosses_boundary = false;
            for j = 1:numel(weight_blocks)
                indices = weight_blocks{j};
                for k = 1:numel(eigen_blocks)
                    overlap = intersect(indices, eigen_blocks{k});
                    if ~isempty(overlap) && numel(overlap) ~= numel(eigen_blocks{k})
                        crosses_boundary = true;
                    end
                end
            end
            if ~crosses_boundary
                target = vectors(:, 1:p);
                solutionfun = @(x) stiefel_distance(x, target, weight_blocks);
                metadata.solution_label = '$d(X_k,\mathcal{X}_\star)/\|X_\star\|_F$';
                metadata.solution_description = ['Relative Frobenius distance, minimized over column signs ', ...
                    'and orthogonal rotations within equal-weight blocks.'];
            else
                solutionfun = @(x) stiefel_subspace_error(x, vectors, eigen_blocks, weight_blocks);
                metadata.solution_label = '$e_{\mathrm{subspace}}(X_k)$';
                metadata.solution_description = ['Invariant spectral-subspace diagnostic, normalized by sqrt(p). ', ...
                    'It vanishes precisely on the optimizer set, including eigenvalue multiplicities ', ...
                    'that cross weight-block boundaries or the rank cutoff; it is not a distance ', ...
                    'to one arbitrarily chosen eigenbasis.'];
            end
            metadata.eigenvalue_tolerance = eigen_tolerance;
            metadata.weight_tolerance = weight_tolerance;

        otherwise
            error('rhzls:UnknownExample', 'Comparison metrics support sphere and stiefel.');
    end
    metricfun = @(x, varargin) evaluate_metrics(x, varargin, problem, solutionfun, optimum, objective_scale);
end

function metrics = evaluate_metrics(x, supplied_cost, problem, solutionfun, optimum, ~)
    if isempty(supplied_cost)
        cost = problem.cost(x);
    else
        cost = supplied_cost{1};
    end
    gap = abs(cost - optimum);
    if isfinite(optimum) && isfield(problem, 'objective_error')
        gap = problem.objective_error(x);
    end
    metrics = struct('objective_error', gap, ...
        'solution_error', solutionfun(x));
end

function [vectors, values, tolerance] = sorted_eigensystem(A)
    [vectors, diagonal] = eig((A + A') / 2);
    [values, order] = sort(diag(diagonal), 'ascend');
    vectors = vectors(:, order);
    % Treat only spectrally indistinguishable eigenvalues as a multiplicity.
    tolerance = 32 * eps(max(1, max(abs(values)))) * max(1, sqrt(size(A, 1)));
end

function blocks = equal_value_blocks(values, tolerance)
    blocks = {};
    first = 1;
    while first <= numel(values)
        last = first;
        while last < numel(values) && abs(values(last + 1) - values(first)) <= tolerance
            last = last + 1;
        end
        blocks{end + 1} = first:last; %#ok<AGROW>
        first = last + 1;
    end
end

function value = sphere_distance(x, minimum_space)
    projection = minimum_space * (minimum_space' * x);
    projection_norm = norm(projection);
    if projection_norm > 0
        nearest = projection / projection_norm;
    else
        nearest = minimum_space(:, 1);
    end
    % Direct subtraction stays informative close to machine precision.
    value = norm(x - nearest);
end

function value = stiefel_distance(x, target, blocks)
    squared_distance = 0;
    for k = 1:numel(blocks)
        indices = blocks{k};
        reference_block = target(:, indices);
        current_block = x(:, indices);
        [left, ~, right] = svd(reference_block' * current_block);
        difference = current_block - reference_block * (left * right');
        squared_distance = squared_distance + norm(difference, 'fro') ^ 2;
    end
    value = sqrt(squared_distance / size(x, 2));
end

function value = stiefel_subspace_error(x, vectors, eigen_blocks, weight_blocks)
    squared_error = 0;
    for j = 1:numel(weight_blocks)
        indices = weight_blocks{j};
        block = x(:, indices);
        allowed_indices = [];
        for k = 1:numel(eigen_blocks)
            multiplicity = numel(intersect(indices, eigen_blocks{k}));
            if multiplicity == 0, continue; end
            eigenspace = vectors(:, eigen_blocks{k});
            singular_values = svd(eigenspace' * block);
            expected = zeros(size(singular_values));
            expected(1:multiplicity) = 1;
            squared_error = squared_error + norm(singular_values - expected) ^ 2;
            allowed_indices = [allowed_indices, eigen_blocks{k}]; %#ok<AGROW>
        end
        allowed = vectors(:, allowed_indices);
        outside = block - allowed * (allowed' * block);
        squared_error = squared_error + norm(outside, 'fro') ^ 2;
    end
    value = sqrt(squared_error / size(x, 2));
end

