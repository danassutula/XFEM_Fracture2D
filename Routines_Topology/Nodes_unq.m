function [id_unq,id_rpt,id_r2u] = Nodes_unq(x,tol)

if nargin == 1 % not recommended
    tol = max(hypot(x(:,1),x(:,2)))*1e-14;
end

tol = tol*tol; 

[n,n_dim] = size(x);
id_r2u = zeros(n,1);

for i = 1:n-1
    if ~id_r2u(i)
        
        q = i+1:n; % test points
        L2 = (x(q,1)-x(i,1)).^2;
        
        for k = 2:n_dim % faster wrt q
            L2 = L2 + (x(q,k)-x(i,k)).^2;
        end
        
        id_r2u(q(L2<tol)) = i;
        
    end
end

id_rpt = find(id_r2u);
id_unq = true(n,1);
id_unq(id_rpt) = 0;
id_unq = find(id_unq);

end