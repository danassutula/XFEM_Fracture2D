function [da,ps] = CrackGrowthSol_GradientBasedOptim_fixedLength(Gs,Hs,da)
%
% Gradient-based active set optimisation considering fixed-length crack tip
% extensions.
%
% This is a proof of concept that energy can be minimised even if only 
% fixed length crack increments are supposed. This relies on the fact that
% if da_inc is sufficeintly small, than Gs and Hs are effectivelly 
% constants (from Taylor's series expansion of potential energy about the
% current crack frcture configuration). Thus, the objective is find da, 
% such that it maximises the associated quadratic function: 
% ps(da) = da'*Gs + 0.5*da'*Hs*da, where da^0(:) = da_inc.
%
% da = crack tip extension lengths
% ps = objective function value (to be maxmimised, e.g. energy dissipated)
%
% Gs = inplane energy release rates
% Hs = rates of energy release rates
%
% Note, for non-competing crack growth (i.e. (Gs_i == Gs_j) == false), the 
% initial solutions is not so important, so long as sum(da) is small or 
% Gs_i >> Hs_ij*da_j for all da such that 0 >= da <= da_inc. Here, we will
% assume that the optimal solution is contained in the trial solution, e.g.
% da^0(:) = da_inc.

if isempty(da) || any(da<0) || all(da==0)
    error('The initial solution for "da" is not valid.')
end

if length(da) == 1
    % trial solution: all tips are extended by da_inc
    da = ones(length(Gs),1)*da;
end

while 1
    
    q = find(da>0); % active set
    
    if length(q) == 1
        break
    end
    
    g = Gs+Hs*da; % real energy dissipation gradient
    d = g(q)-mean(g(q)); % advance dir. (instantanious, ideal)
    
    % Assume a crack does not grow at all if it is less energetically 
    % favourable. Thus, the modified crack extensions are updated to:
    
    da(q(d<0)) = 0;
    
end

if sum(da) == 0
   error('Zero solution?')
end

% the resulting value of the objective function
ps = da'*Gs+0.5*(da'*Hs*da);

end
