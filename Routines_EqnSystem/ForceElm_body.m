function Fe = ForceElm_body(x_elm,N,dNdX,W,f)

%--------------------------------------------------------------------------

% Fe = int N*f dOmg

% (!!!) if the element is enriched, the enrichment must remain constant

%--------------------------------------------------------------------------

nJ = size(x_elm,1);
nN = size(N,2);

Fe = zeros(2*nN,1);
S  = zeros(1,nN);

qj = 1:nJ;
qd = 1:2;

for i = 1:length(W)

    J = dNdX(qd,qj)*x_elm;
    A = (J(1)*J(4)-J(2)*J(3))*W(i);
    
    S  = S  + A*N(i,:);
    qd = qd + 2;

end

Fe(1:2:end-1) = S*f(1);
Fe(2:2:end)   = S*f(2);
    
end
