function matrix = full_matrix(x)
%FULL_MATRIX Materialize a numeric array or a Manopt U,S,V triplet.
%   Intended for small examples and tests, not large-scale optimization.
if isnumeric(x)
    matrix = x;
elseif isstruct(x) && all(isfield(x, {'U', 'S', 'V'}))
    matrix = x.U*x.S*x.V';
elseif isstruct(x) && all(isfield(x, {'G', 'H'}))
    matrix = x.G*x.H';
else
    error('rhzls:InvalidRepresentation', 'Expected a numeric matrix, U,S,V triplet, or G,H factors.');
end
end
