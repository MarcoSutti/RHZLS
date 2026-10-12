function [stepsize, newx, newkey, lsstats] = manopt_linesearch( ...
        problem, x, d, f0, df0, options, storedb, key)
%MANOPT_LINESEARCH Adapter for options.linesearch in Manopt solvers.
% options.rhzls is the RHZLS parameter struct (or problem.options_hzls).
% Returns stepsize = alpha*norm(d), as required by Manopt.
user = struct();
if isfield(problem, 'options_hzls'), user = problem.options_hzls; end
if isfield(options, 'rhzls'), user = options.rhzls; end
user = rhzls.options(user);
if isfield(storedb.internal, 'rhzls_alpha')
    user.alpha0 = min(2*storedb.internal.rhzls_alpha, user.alphamax);
end
[alpha, newx, lsstats] = rhzls.linesearch(problem, x, d, user, f0, df0);
stepsize = alpha*problem.M.norm(x, d);
lsstats.alpha = alpha;
lsstats.stepsize = stepsize;
if lsstats.accepted
    storedb.internal.rhzls_alpha = alpha;
    newkey = storedb.getNewKey();
else
    error('rhzls:LineSearchFailed', ...
        'RHZLS failed (%s). No Wolfe step was accepted.', lsstats.status);
end
end
