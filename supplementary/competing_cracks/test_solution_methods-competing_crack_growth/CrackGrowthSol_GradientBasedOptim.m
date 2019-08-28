function [da,ps] = CrackGrowthSol_GradientBasedOptim(Gs,Hs,da)
%
% ActiveSetOptim_gradient for the off-line solution of competing crack
% growth where the potential energy (Pi) of the solid is approximated by a
% quadratic: $Pi(Da) = Pi0 - Gs_i*Da_i - 1/2 Hs_ij*Da_i*Da_j$
%
% maximize ps(da) = (da'*Gs + 1/2*da'*Hs*da),
% subject to sum(da) = da_tot and da >= 0.
%
% NOTE: the difference between "_discrete" and non-discrete is that in
% the discrete form, cracks are competing over the finite increment in 
% the total fracture area 'da_tot'. On the other hand, in the non-discrete
% form of the active set gradient-based optimisation, the cracks are 
% assumed to be competting prior to the crack extension 'da_tot'; in other
% words, 'Gs' is equal at all the competing crack tips and, consequently,
% Gs can be ignored as all the essential information is in 'Hs'
%
% An energy minimization algorithm based on a gradient-descent scheme 
% subject to the constraint of a fixed total fracture increment $da_tot$
%
% da = fracture increment values / initial trial solution
% ps = objective function value
% gs = mean energy release rate
%
% Gs = inplane energy release rates
% Hs = rates of energy release rates
% da_tot = L1-norm of fracture extension
%
% n.b. $Hs$ will be scaled by $da_tot$ so that $da$ is relative, i.e.
% norm(da,1) = 1; however, in the end, $da$ is scaled back to $da_tot$.

da_tot = sum(da);
da = da/da_tot;

n = size(Hs,1);
e = ones(n,1);

% finding reduced Hs
Hs = Hs.*da_tot; % (scaling Hs)
C = eye(n)-(e*e')/n;
Hs_red = C*Hs*C;

% reduced eig. val.
u_red = eig(Hs_red);
u_inf = norm(u_red,inf);

tol = 1e-6;

if all(u_red < tol*u_inf)
    
%     fprintf('\nReduced-Hessian is negative semidefinite.\n');
    
elseif all(u_red > -tol*u_inf)
    
%     fprintf('\nReduced-Hessian is positive semidefinite.\n');
    % -> only one crack can grow (assuming \Delta a_i = 1 or 0)
    
    [ps,i] = max(Gs+0.5*diag(Hs));
    
    da = zeros(n,1);
    da(i) = da_tot;
    ps = ps*da_tot;
    
    return
    
end

% current energy value
ps = da'*Gs+0.5*(da'*Hs*da);

while 1
    
    da0 = da;
    ps0 = ps;
    
    d = zeros(n,1);
    p = true(n,1);
    g = Gs+Hs*da;
    
    while 1
        
        % direction satisfies L_1-norm constraint
        % gd(p) = g(p)-(e(p)*e(p)')/(e(p)'*e(p))*g(p);
        
        d(p) = g(p)-mean(g(p)); % advance direction
        q = da < tol & d < 0; % assess feasibility
        
        % since only the smallest members can be removed each time, the 
        % mean value has to increase. As a result, any previously negative 
        % points will remain negative. For efficiency, remove all negative 
        % points at once (i.e. points that are below the mean).
        
        if any(q)
            d(q) = 0;
            p(q) = 0;
        else
            break
        end
        
    end
    
    % could be infeasibile
    q = da > 0 & d < 0;
    
    if ~any(q)
        break
    end
    
    % enforce feasibility: ai+w*di >= 0 
    w = min(da(q)./-d(q)); 
        
    if d'*Hs*d < 0 % (concave down)  
        % try line-search for an optimum w   
        w = min(w,-(d'*Gs+d'*Hs*da)/(d'*Hs*d));
    end
    
    da = da + d*w;
    
    if abs(sum(da,1)-1) > tol
        % precision lost
        if any(da < tol)
            da = da0; % failed: use prev. solution then stop
        else % re-normalize (L1); attempt to continue anyhow
            da = da/norm(da,1);
        end
    end
    
    ps = da'*Gs+0.5*(da'*Hs*da);
    
    if ps <= ps0
        ps=ps0; da=da0;
        break % converged
    end
    
end

da = da*da_tot;
ps = ps*da_tot;

end
