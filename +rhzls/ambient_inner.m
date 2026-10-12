function value = ambient_inner(gradient, velocity)
%AMBIENT_INNER Real Frobenius pairing of ambient gradient and curve velocity.
%   Numeric matrices (including sparse ones) and U,S,V triplets are
%   supported for GRADIENT. VELOCITY may be numeric or factored as G,H,
%   meaning G*H'. The factored path never builds a dense ambient matrix.
if isnumeric(velocity)
    if isnumeric(gradient)
        if ~isequal(size(gradient), size(velocity))
            error('rhzls:InvalidGradient', 'Ambient gradient and velocity sizes must agree.');
        end
        value = real(gradient(:)'*velocity(:));
    elseif isstruct(gradient) && all(isfield(gradient, {'U', 'S', 'V'}))
        value = real(trace(gradient.S'*(gradient.U'*velocity*gradient.V)));
    else
        error('rhzls:InvalidGradient', 'Unsupported ambient gradient representation.');
    end
elseif isstruct(velocity) && all(isfield(velocity, {'G', 'H'}))
    if isnumeric(gradient)
        product = gradient*velocity.H;
        value = real(sum(sum(conj(velocity.G).*product)));
    elseif isstruct(gradient) && all(isfield(gradient, {'U', 'S', 'V'}))
        value = real(trace((velocity.G'*gradient.U)*gradient.S*(gradient.V'*velocity.H)));
    else
        error('rhzls:InvalidGradient', 'Unsupported ambient gradient representation.');
    end
else
    error('rhzls:InvalidGradient', 'Velocity must be numeric or have G,H factor fields.');
end
end
