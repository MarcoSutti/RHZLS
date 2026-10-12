function [M, curve, maxstep] = manifold(name, dimensions)
%MANIFOLD Manopt geometry with an explicitly matched retraction velocity.
%   [M, CURVE, MAXSTEP] = RHZLS.MANIFOLD(NAME, DIMENSIONS) supports
%   'sphere' with n, 'stiefel' with [n p], and 'fixedrank' with [m n r].
%   CURVE(x, d, t) returns the retracted point and its ambient velocity.
%   MAXSTEP(x, d) bounds the positive smooth part of the retraction curve.

if isstring(name) && isscalar(name)
    name = char(name);
end
if ~ischar(name) || ~isvector(dimensions) || isempty(dimensions) || ...
        any(~isfinite(dimensions)) || any(dimensions < 1) || ...
        any(dimensions ~= round(dimensions))
    error('rhzls:InvalidGeometry', 'Supply a manifold name and positive integer dimensions.');
end
dimensions = dimensions(:).';
name = lower(strtrim(name));
switch name
    case 'sphere'
        require_factory('spherefactory');
        if numel(dimensions) ~= 1 || dimensions(1) < 2
            error('rhzls:InvalidGeometry', 'Sphere dimensions must be a scalar n >= 2.');
        end
        M = spherefactory(dimensions(1));
        curve = @rhzls.sphere_curve;
        maxstep = @(x, d) Inf;
    case 'stiefel'
        require_factory('stiefelfactory');
        if numel(dimensions) ~= 2 || dimensions(1) < dimensions(2)
            error('rhzls:InvalidGeometry', 'Stiefel dimensions must be [n p], with n >= p.');
        end
        M = stiefelfactory(dimensions(1), dimensions(2));
        curve = @rhzls.stiefel_curve;
        maxstep = @(x, d) Inf;
        if isfield(M, 'invretr_polar')
            M.invretr = M.invretr_polar;
        end
    case 'fixedrank'
        require_factory('fixedrankembeddedfactory');
        if numel(dimensions) ~= 3 || dimensions(3) > min(dimensions(1:2))
            error('rhzls:InvalidGeometry', 'Fixed-rank dimensions must be [m n r], with r <= min(m,n).');
        end
        % RMGLS bundled a four-argument modified Manopt factory. Stock
        % Manopt needs only three arguments. In both cases we replace its
        % retraction with the same exact AO curve differentiated below.
        if nargin('fixedrankembeddedfactory') == 4
            M = fixedrankembeddedfactory(dimensions(1), dimensions(2), ...
                                         dimensions(3), 'ortho_ao');
        else
            M = fixedrankembeddedfactory(dimensions(1), dimensions(2), dimensions(3));
        end
        curve = @rhzls.fixedrank_curve;
        maxstep = @rhzls.fixedrank_maxstep;
    otherwise
        error('rhzls:InvalidGeometry', 'Unknown manifold "%s". Use sphere, stiefel, or fixedrank.', name);
end
M.retr = curve;
M.rhzls_maxstep = maxstep;
end

function require_factory(name)
if exist(name, 'file') ~= 2
    error('rhzls:MissingManopt', 'Manopt factory %s is missing. Run startup_rhzls first.', name);
end
end
