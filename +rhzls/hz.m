function [alpha, f, df, stats] = hz(oracle, user, f0, df0)
%HZ Scalar Hager-Zhang search: bracketing, secant^2, and bisection.
% oracle(t) returns [phi(t), phi'(t)].
% Accepts ordinary OR approximate Wolfe conditions.
% A failure returns alpha=0, the origin values, and stats.accepted=false.
% No Armijo-only fallback is silently accepted.
% Reference: Hager & Zhang, Algorithm 851, ACM TOMS 32 (2006), section 3.
if nargin < 2, user = struct(); end
o = rhzls.options(user);
alpha = 0;
stats = struct('accepted', false, 'status', 'max_evaluations', ...
    'evaluations', 0, 'iterations', 0, 'wolfe', false, ...
    'approxwolfe', false, 'phi_limit', NaN, 'bracket', [NaN NaN], ...
    'f', NaN, 'df', NaN, 'f0', NaN, 'df0', NaN);
if nargin < 4
    [f0, df0] = oracle(0);
    stats.evaluations = 1;
end
validate_value(f0); validate_value(df0);
f = f0; df = df0;
stats.f0 = f0; stats.df0 = df0;
stats.f = f0; stats.df = df0;
if ~isfinite(f0) || ~isfinite(df0)
    stats.status = 'invalid_origin'; return;
end
if df0 >= 0
    stats.status = 'not_descent'; return;
end
limit = f0 + o.epsilon*max(abs(f0), o.function_scale);
stats.phi_limit = limit;
done = false;
best = point(0, f0, df0);
a = best; b = best;
c = min(o.alpha0, o.alphamax);

% B0-B3: expand until the derivative changes sign. A high objective with
% negative slope invokes U3, rather than accepting a spurious later valley.
while ~done
    p = trial(c);
    if done, break; end
    if finite(p) && p.df >= 0
        b = p; break;
    elseif ~finite(p) || p.f > limit
        [a, b] = bisect_high(a, p);
        break;
    else
        a = p;
        if c >= o.alphamax
            stop('alphamax'); break;
        end
        next = min(o.rho*c, o.alphamax);
        if next <= c
            stop('interval_too_small'); break;
        end
        c = next;
    end
end

% L1-L3. update preserves phi(a)<=limit, phi'(a)<0<=phi'(b).
while ~done
    stats.bracket = [a.t b.t];
    if ~(b.t > a.t && finite(a) && finite(b) && a.df < 0 && b.df >= 0)
        stop('invalid_bracket'); break;
    end
    if too_small(a.t, b.t)
        stop('interval_too_small'); break;
    end
    stats.iterations = stats.iterations + 1;
    olda = a; oldb = b;
    c = secant(a, b);
    if ~interior(c, a.t, b.t), c = midpoint(a.t, b.t); end
    p = trial(c);
    if done, break; end
    [a, b] = update(a, b, p);
    if done, break; end

    % S2-S3: interpolate using the changed endpoint and its OLD value.
    % This is different from repeating a secant on the new full interval.
    c2 = NaN;
    if b.t == p.t
        c2 = secant(oldb, b);
    elseif a.t == p.t
        c2 = secant(olda, a);
    end
    if interior(c2, a.t, b.t)
        p = trial(c2);
        if done, break; end
        [a, b] = update(a, b, p);
        if done, break; end
    end
    if b.t-a.t > o.gamma*(oldb.t-olda.t)
        p = trial(midpoint(a.t, b.t));
        if done, break; end
        [a, b] = update(a, b, p);
    end
end
if stats.accepted
    alpha = best.t; f = best.f; df = best.df;
end
stats.f = f; stats.df = df;

    function p = trial(t)
        if stats.evaluations >= o.maxevals
            stop('max_evaluations'); p = point(t, NaN, NaN); return;
        end
        [ft, dft] = oracle(t);
        validate_value(ft); validate_value(dft);
        stats.evaluations = stats.evaluations + 1;
        p = point(t, ft, dft);
        if finite(p) && t > 0
            w = ft <= f0 + o.delta*t*df0 && dft >= o.sigma*df0;
            aw = ft <= limit && dft >= o.sigma*df0 ...
                && dft <= (2*o.delta-1)*df0;
            if (strcmp(o.acceptance, 'either') && w) || aw
                best = p;
                stats.accepted = true;
                stats.wolfe = w; stats.approxwolfe = aw;
                if w && strcmp(o.acceptance, 'either'), stats.status = 'wolfe'; else, stats.status = 'approximate_wolfe'; end
                done = true;
            end
        end
    end

    function [lo, hi] = update(lo, hi, p)
        if ~finite(p)
            [lo, hi] = bisect_high(lo, p);
        elseif p.df >= 0
            hi = p;
        elseif p.f <= limit
            lo = p;
        else
            [lo, hi] = bisect_high(lo, p);
        end
    end

    function [lo, hi] = bisect_high(lo, hi)
        % U3: high value with negative slope, or an invalid domain trial.
        while ~done
            if too_small(lo.t, hi.t)
                stop('interval_too_small'); return;
            end
            t = lo.t + o.theta*(hi.t-lo.t);
            p = trial(t);
            if done, return; end
            if finite(p) && p.df >= 0
                hi = p; return;
            elseif finite(p) && p.f <= limit
                lo = p;
            else
                hi = p;
            end
        end
    end

    function stop(status)
        stats.status = status; done = true;
    end
end

function p = point(t, f, df)
p = struct('t', t, 'f', f, 'df', df);
end
function yes = finite(p)
yes = isfinite(p.f) && isfinite(p.df);
end
function yes = interior(t, a, b)
yes = isfinite(t) && t > a && t < b;
end
function yes = too_small(a, b)
yes = b-a <= 4*eps(max(abs([a b])));
end
function c = midpoint(a, b)
c = a + .5*(b-a);
end
function c = secant(a, b)
denom = b.df-a.df;
if denom == 0 || ~isfinite(denom)
    c = NaN;
else
    c = a.t - a.df*((b.t-a.t)/denom);
end
end
function validate_value(v)
if ~isnumeric(v) || ~isreal(v) || ~isscalar(v)
    error('rhzls:Oracle', 'The oracle must return real scalar cost and derivative.');
end
end
