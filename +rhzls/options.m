function o = options(user)
%OPTIONS Validated parameters for the scalar Hager-Zhang line search.
% phi_limit = phi(0) + epsilon*max(abs(phi(0)), function_scale).
% maxevals bounds ALL oracle calls, including an uncached evaluation at zero.
% acceptance='either' (default) permits Wolfe or approximate Wolfe;
% acceptance='approximate' requires both approximate-Wolfe derivative bounds.
o = struct('delta', .1, 'sigma', .9, 'epsilon', 1e-6, ...
    'function_scale', 1, 'alpha0', 1, 'alphamax', 1e6, ...
    'rho', 5, 'gamma', .66, 'theta', .5, 'maxevals', 100, 'acceptance', 'either');
if nargin > 0 && ~isempty(user)
    if ~isstruct(user) || ~isscalar(user)
        error('rhzls:Options', 'Options must be a scalar struct.');
    end
    names = fieldnames(user);
    for k = 1:numel(names)
        if ~isfield(o, names{k})
            error('rhzls:Options', 'Unknown RHZLS option: %s.', names{k});
        end
        o.(names{k}) = user.(names{k});
    end
end
if ~(ischar(o.acceptance) || (isstring(o.acceptance) && isscalar(o.acceptance)))
    error('rhzls:Options', 'acceptance must be either or approximate.');
end
o.acceptance = char(o.acceptance);
if ~any(strcmp(o.acceptance, {'either','approximate'}))
    error('rhzls:Options', 'acceptance must be either or approximate.');
end
names = setdiff(fieldnames(o), {'acceptance'});
for k = 1:numel(names)
    v = o.(names{k});
    if ~isnumeric(v) || ~isreal(v) || ~isscalar(v) || ~isfinite(v)
        error('rhzls:Options', '%s must be a finite real scalar.', names{k});
    end
end
if ~(o.delta > 0 && o.delta < .5 && o.sigma >= o.delta && o.sigma < 1 ...
        && o.epsilon >= 0 && o.function_scale >= 0 ...
        && o.alpha0 > 0 && o.alphamax > 0 && o.rho > 1 ...
        && o.gamma > 0 && o.gamma < 1 && o.theta > 0 && o.theta < 1 ...
        && o.maxevals >= 1 && o.maxevals == floor(o.maxevals))
    error('rhzls:Options', 'Invalid RHZLS parameter range. See help rhzls.options.');
end
end
