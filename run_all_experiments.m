function results = run_all_experiments(show_figures)
%RUN_ALL_EXPERIMENTS Repeat all benchmarks, save results, and export nine PDFs.
%   results = run_all_experiments() displays the nine final figures.
%   results = run_all_experiments(false) exports PDFs without opening windows.
%   The caller's current folder and defaultFigureVisible setting are restored.
%   Results and PDFs are written to this package's docs directory.
%
%   The returned structure contains demo, accuracy, and candidates results,
%   output paths, and aggregate line-search function/gradient counts.

    if nargin < 1 || isempty(show_figures), show_figures = true; end
    validateattributes(show_figures, {'logical', 'numeric'}, ...
        {'scalar', 'real', 'finite'}, mfilename, 'show_figures');
    if ~ismember(show_figures, [0, 1])
        error('rhzls:InvalidOptions', 'show_figures must be true or false.');
    end

    project_root = fileparts(mfilename('fullpath'));
    caller_folder = pwd;
    previous_visibility = char(get(groot, 'defaultFigureVisible'));
    restore_caller = environment_cleanup(caller_folder, previous_visibility); %#ok<NASGU>
    cd(project_root);
    paths = startup_rhzls();
    paths.docs = fullfile(project_root, 'docs');
    if exist(paths.docs, 'dir') ~= 7, mkdir(paths.docs); end
    set(groot, 'defaultFigureVisible', 'off');

    fprintf('\nRHZLS: running all benchmark suites.\n');
    [x, f, info, details] = run_rhzls('all', struct('plot', false));
    save(fullfile(paths.docs, 'demo_results.mat'), 'x', 'f', 'info', 'details', 'paths');
    demo = struct('x', x, 'cost', f, 'info', info, 'details', details);

    accuracy = run_accuracy_examples(false);
    require_success(accuracy, {'manopt', 'rhzls', 'rhzls_approx'}, 'Accuracy');
    candidates = explore_manopt_examples(false);
    require_success(candidates, {'manopt', 'rhzls'}, 'Manopt-derived');

    % Export the complete set once with the requested visibility, so no
    % intermediate or duplicate figure windows appear.
    if show_figures, visibility = 'on'; else, visibility = 'off'; end
    set(groot, 'defaultFigureVisible', visibility);
    regenerate_plots(fullfile(paths.docs, 'plot_validation.m'));
    pdf_names = {'comparison_stiefel.pdf', 'comparison_sphere.pdf', ...
        'accuracy_sphere_positive.pdf', 'accuracy_stiefel_positive.pdf', ...
        'accuracy_orthogonal_procrustes.pdf', 'accuracy_rank_approximation.pdf', ...
        'candidate_generalized_eigen.pdf', 'candidate_thomson_tetrahedron.pdf', ...
        'candidate_maxcut_cycle5.pdf'};
    paths.pdfs = cellfun(@(name) fullfile(paths.docs, name), pdf_names, ...
        'UniformOutput', false);
    for k = 1:numel(paths.pdfs)
        if exist(paths.pdfs{k}, 'file') ~= 2
            error('rhzls:MissingOutput', 'Expected PDF was not generated: %s', paths.pdfs{k});
        end
    end
    paths.mat_files = {fullfile(paths.docs, 'demo_results.mat'), ...
        fullfile(paths.docs, 'accuracy_examples.mat'), ...
        fullfile(paths.docs, 'manopt_candidate_results.mat')};

    totals = count_evaluations(demo, accuracy, candidates);
    results = struct('demo', demo, 'accuracy', accuracy, ...
        'candidates', candidates, 'paths', paths, 'totals', totals);
    fprintf('\nCompleted %d solver runs; saved three MAT files and nine PDFs in:\n  %s\n', ...
        totals.runs, paths.docs);
    fprintf('Line-search totals: nf = %d, ng = %d (all runs, including accuracy controls).\n', ...
        totals.nf, totals.ng);
    % Restore explicitly before returning function handles in the results.
    % Keep exception cleanup in its own workspace to avoid delayed lifetime.
    restore_environment(caller_folder, previous_visibility);
    clear restore_caller;
end

function require_success(entries, methods, suite)
    if isempty(entries)
        error('rhzls:ExperimentFailed', '%s suite returned no experiments.', suite);
    end
    failures = {};
    for k = 1:numel(entries)
        for j = 1:numel(methods)
            method = methods{j};
            if ~isfield(entries(k), method)
                reason = 'missing method result';
            elseif ~isempty(entries(k).(method).error)
                reason = entries(k).(method).error;
            elseif isempty(entries(k).(method).info)
                reason = 'empty iteration history';
            else
                continue;
            end
            failures{end+1} = sprintf('%s, seed %d, %s: %s', ...
                entries(k).family, entries(k).seed, method, reason); %#ok<AGROW>
        end
    end
    if ~isempty(failures)
        error('rhzls:ExperimentFailed', ...
            '%s suite failed. Review saved results before using its PDFs:\n%s', ...
            suite, strjoin(failures, newline));
    end
end

function totals = count_evaluations(demo, accuracy, candidates)
    histories = {};
    for example = {'stiefel', 'sphere'}
        comparison = demo.details.(example{1}).comparison;
        histories{end+1} = comparison.manopt.info; %#ok<AGROW>
        histories{end+1} = comparison.rhzls.info; %#ok<AGROW>
    end
    suites = {accuracy, candidates};
    methods = {{'manopt', 'rhzls', 'rhzls_approx'}, {'manopt', 'rhzls'}};
    for suite = 1:numel(suites)
        entries = suites{suite};
        for k = 1:numel(entries)
            for j = 1:numel(methods{suite})
                histories{end+1} = entries(k).(methods{suite}{j}).info; %#ok<AGROW>
            end
        end
    end
    totals = struct('runs', numel(histories), 'nf', 0, 'ng', 0);
    for k = 1:numel(histories)
        history = histories{k};
        totals.nf = totals.nf + sum([history.nf]);
        totals.ng = totals.ng + sum([history.ng]);
    end
end

function regenerate_plots(script_path)
    % Scripts use their own variables; keep them out of the runner workspace.
    run(script_path);
end

function restore_environment(folder, visibility)
    set(groot, 'defaultFigureVisible', visibility);
    cd(folder);
end

function cleanup = environment_cleanup(folder, visibility)
    cleanup = onCleanup(@() restore_environment(folder, visibility));
end
