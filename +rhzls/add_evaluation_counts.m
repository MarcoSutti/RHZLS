function info = add_evaluation_counts(info, line_search)
%ADD_EVALUATION_COUNTS Actual line-search callback counts for each iteration.
% nf(k), ng(k) exclude outer-solver evaluations and diagnostics. Iteration
% zero has no line search, so both counts are zero there. Oracle attempts
% are a different quantity and cannot substitute for callback counts.
% Official Armijo records costevals but no gradevals. For our separate-cost
% benchmark contract it makes no gradient calls; permit ng=0 only when the
% actual solver-selected callback is the official @linesearch.
if nargin < 2, line_search = []; end
for k = 1:numel(info)
    info(k).nf = 0;
    info(k).ng = 0;
    if info(k).iter == 0, continue; end
    if ~isfield(info(k), 'linesearch') ...
            || ~isfield(info(k).linesearch, 'costevals')
        error('rhzls:MissingEvaluationCounts', ...
            ['Saved histories lack actual evaluation counts. ' ...
             'Run run_all_experiments(false) to regenerate them.']);
    end
    nf = info(k).linesearch.costevals;
    if isfield(info(k).linesearch, 'gradevals')
        ng = info(k).linesearch.gradevals;
    elseif isequal(line_search, @linesearch)
        ng = 0;
    else
        error('rhzls:MissingEvaluationCounts', ...
            'Gradient counts are absent and the callback is not official Armijo.');
    end
    validateattributes(nf, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
    validateattributes(ng, {'numeric'}, {'scalar', 'real', 'finite', 'integer', 'nonnegative'});
    info(k).nf = nf;
    info(k).ng = ng;
end
end
