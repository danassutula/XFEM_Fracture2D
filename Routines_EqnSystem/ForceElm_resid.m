function Fe = ForceElm_resid(x_elm,dNdX,W,f)

%--------------------------------------------------------------------------

% Fe = - int B*(D*eps_0 + sig_0) dOmg = - int B*(f) dOmg

% (!!!) if the element is enriched, the enrichment must remain constant

%--------------------------------------------------------------------------

f = -f; % (since RHS is negative)

nJ = size(x_elm,1);
nN = size(dNdX,2);

Fe = zeros(2*nN,1);
B  = zeros(2,nN);
C  = zeros(2,2);

qj = 1:nJ;
qd = 1:2;

for w = W(:)'
    
    J = dNdX(qd,qj)*x_elm;
    
    % C = inv(J)*det(J) = adj(J)
    C(1) =  J(4); C(3) = -J(3); 
    C(2) = -J(2); C(4) =  J(1);
    
    B = B  + (C*dNdX(qd,:)).*w;
    
    qd = qd + 2;
    
end

Fe(1:2:end-1) = B(1,:)*f(1)+B(2,:)*f(3);
Fe(2:2:end)   = B(2,:)*f(2)+B(1,:)*f(3);

end
