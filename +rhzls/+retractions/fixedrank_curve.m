function [y, velocity] = fixedrank_curve(x, d, t)
%FIXEDRANK_CURVE Absil-Oseledets orthographic retraction and exact velocity.
%   Points use Manopt's U,S,V representation and tangents use Up,M,Vp.
%   For K=S+t*M, L=U*K+t*Up and R=V*K'+t*Vp, the curve is L*(K\R').
%   The velocity is returned as factors G,H with velocity=G*H'. Neither
%   the retraction nor the derivative constructs an m-by-n dense matrix.
%   This local retraction requires K to be nonsingular. Singular trials
%   are rejected instead of modifying K or its returned singular values.
if nargin < 3
    t = 1;
end
if ~isscalar(t) || ~isreal(t) || ~isfinite(t) || ...
        ~isstruct(x) || ~all(isfield(x, {'U', 'S', 'V'})) || ...
        ~isstruct(d) || ~all(isfield(d, {'Up', 'M', 'Vp'}))
    error('rhzls:InvalidRetraction', 'Fixed-rank retraction requires a point and a Manopt tangent.');
end
r = size(x.S, 1);
if r == 0 || ~isequal(size(x.S), [r r]) || ~isequal(size(d.M), [r r]) || ...
        size(x.U, 2) ~= r || size(x.V, 2) ~= r || ...
        ~isequal(size(d.Up), size(x.U)) || ~isequal(size(d.Vp), size(x.V))
    error('rhzls:InvalidRetraction', 'Incompatible fixed-rank point and tangent dimensions.');
end
values = {x.U, x.S, x.V, d.Up, d.M, d.Vp};
for j = 1:numel(values)
    if ~isnumeric(values{j}) || ~isreal(values{j}) || any(~isfinite(values{j}(:)))
        error('rhzls:InvalidRetraction', 'Fixed-rank inputs must be finite real matrices.');
    end
end
k = x.S + t*d.M;
if any(~isfinite(k(:)))
    error('rhzls:InvalidRetraction', 'Nonfinite orthographic retraction core.');
end
if rcond(k) <= 1e-12
    error('rhzls:SingularRetraction', 'The orthographic retraction core S+t*M is numerically singular.');
end
if t == 0 || (all(d.M(:) == 0) && all(d.Up(:) == 0) && all(d.Vp(:) == 0))
    % Preserve R_x(0)=x exactly, including its factor convention.
    y = x;
    if nargout > 1
        velocity.G = [x.U*d.M + d.Up, x.U];
        velocity.H = [x.V, d.Vp];
    end
    return;
end
left = x.U*k + t*d.Up;
right = x.V*k.' + t*d.Vp;
[ql, rl] = qr(left, 0);
[qr_, rr] = qr(right, 0);
core = rl*(k\rr.');
if any(~isfinite(core(:)))
    error('rhzls:InvalidRetraction', 'Nonfinite orthographic retraction trial.');
end
[u, s, v] = svd(core, 'econ');
if any(diag(s) <= 0) || any(~isfinite(s(:)))
    error('rhzls:InvalidRetraction', 'The orthographic trial does not have the prescribed rank.');
end
y.U = ql*u;
y.S = s;
y.V = qr_*v;
if nargout > 1
    f = left/k;
    h = k\right.';
    leftdot = x.U*d.M + d.Up;
    rightdot = x.V*d.M.' + d.Vp;
    velocity.G = [leftdot, f, -f*d.M];
    velocity.H = [h.', rightdot, h.'];
    if any(~isfinite(velocity.G(:))) || any(~isfinite(velocity.H(:)))
        error('rhzls:InvalidRetraction', 'Nonfinite orthographic retraction velocity.');
    end
end
end
