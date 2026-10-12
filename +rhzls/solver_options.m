function o = solver_options(user)
%SOLVER_OPTIONS Shared stopping criteria for standalone SD and comparisons.
% Zero disables the gradient stopping test for plateau experiments.
o = struct('maxiter', 3000, 'tolgradnorm', 0, ...
    'minstepsize', 0, 'verbosity', 1, 'linesearch', struct());
if nargin < 1 || isempty(user), user = struct(); end
if ~isstruct(user) || ~isscalar(user)
    error('rhzls:Options', 'Solver options must be a scalar struct.');
end
names = fieldnames(user);
for k = 1:numel(names)
    if ~isfield(o, names{k})
        error('rhzls:Options', 'Unknown solver option: %s.', names{k});
    end
    o.(names{k}) = user.(names{k});
end
validateattributes(o.maxiter, {'numeric'}, {'scalar', 'integer', 'nonnegative', 'finite'});
validateattributes(o.tolgradnorm, {'numeric'}, {'scalar', 'nonnegative', 'finite'});
validateattributes(o.minstepsize, {'numeric'}, {'scalar', 'nonnegative', 'finite'});
validateattributes(o.verbosity, {'numeric'}, {'scalar', 'nonnegative', 'finite'});
rhzls.options(o.linesearch);
end
