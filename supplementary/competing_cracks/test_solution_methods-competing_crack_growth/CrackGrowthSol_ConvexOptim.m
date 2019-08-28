function [da,ps] = CrackGrowthSol_ConvexOptim(Gs,Hs,da_tot)

% use for a convex (up or down) objective function: 
% ps(da) = da'*(Gs+1/2*Hs*da), i.e. Hs_red = C*Hs*C
% needs to be pos/neg semidefinite

% da = crack tip increment values
% ps = objective function value (to be maximised)

% Gs = inplane energy release rates
% Hs = rates of energy release rates

% da_tot = total fracture extension
% da_tot == sum(da)

% (!) Here, Hs is scaled by da_tot so that da is relative, i.e.
% sum(da) = 1; however, in the end, da is scaled back to da_tot.

tol = 1e-9;

n = size(Hs,1);
e = ones(n,1);

% determine the behaviour of Hs in the reduced space of sum(da) == da_tot
Hs = Hs.*da_tot;
C = eye(n)-(e*e')/n;
Hs_red = C*Hs*C;

% reduced eig. val.
u_red = eig(Hs_red);
u_inf = norm(u_red,inf);

if all(u_red < tol*u_inf)
    
    fprintf('\nReduced-Hessian is negative semidefinite.\n')
    
elseif all(u_red > -tol*u_inf)
    
    fprintf('\nReduced-Hessian is positive semidefinite.\n')
    % -> only one crack can grow (assuming \Delta a_i = 1 or 0)
    
    [ps,i] = max(Gs+0.5*diag(Hs));
    
    da = zeros(n,1);
    da(i) = da_tot;
    ps = ps*da_tot;
    
    return
    
else
    warning('Reduced-Hessian is indefinite. The solution could be wrong.')
end

% gradient of residual
dR = [Hs,e;e',0];

n_trials = 1;

Da = zeros(n,n_trials);
Ps = zeros(1,n_trials);

for i = 1:n_trials

lm = 0; % Lg. multiplier

if n_trials > 1
    % initial solution
    da = ones(n,1);
    da = da/sum(da);
else
    da = abs(rand(n,1));
    da = da/sum(da);
end
    
ps = da'*Gs+0.5*(da'*Hs*da); 

% remember this:
% Hs <- Hs*da_tot,
% ps <- ps/da_tot

while 1
    
    p = true(n+1,1); % +1 for the lagrange multiplier
    d = zeros(n+1,1); % +1 for the lagrange multiplier
    
    g = Gs+Hs*da+e*lm;
    c = e'*da-1;
    
    % residual
    R = [g;c];
    
    while 1
        
        d(p) = -dR(p,p)\R(p);
        q = [da<tol;false] & d<0;
        
        if nnz(q) > 1 % SAFETY MEASURE
            g = Gs+Hs*da; r = find(q);
            q(r(g(r)>min(g(r)))) = 0;
            % (do not discard multiple)
        end
        
        if any(q)
            d(q) = 0;
            p(q) = 0;
        else
            break
        end
        
    end
    
    if norm(d(1:end-1)) < tol
        break
    end
    
    % could be infeasibile
    q = [da>0;false] & d<0;
    
    if any(q) % enforce feasibility: ai+w*di >= 0
        w = min([1;da(q(1:end-1))./-d(q)]);
    else 
        w = 1;
    end
    
    lm = lm + d(end)*w;
    da = da + d(1:end-1)*w;
    ps = da'*Gs+0.5*(da'*Hs*da);
    
end

Da(:,i) = da;
Ps(1,i) = ps;

end

[ps,i] = max(Ps);
da = Da(:,i)*da_tot;
ps = ps*da_tot;

end
