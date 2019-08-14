
%--------------------------------------------------------------------------
% Constitutive data
%--------------------------------------------------------------------------

if nPhase > length(E) || nPhase > length(v)
    error('Undefined material phases, n_phz = %d',nPhase)
end

% keep in column 
E = E(1:nPhase); E = E(:);
v = v(1:nPhase); v = v(:);

% constitutive matrix (Hooke's law)
[cDMatx,alpha,kappa] = MtxConstit(nPhase,E,v,problemType);

% for SIF calc.
E_str = E./alpha;

% shear modulus
G = E./(1+v)/2;

%--------------------------------------------------------------------------