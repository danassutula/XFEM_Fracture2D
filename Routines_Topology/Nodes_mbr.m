function [id_mbr,id_ref] = Nodes_mbr(x_mbr,x_ref,tol)
%--------------------------------------------------------------------------

% tol = eps('double')*uElSiz

%--------------------------------------------------------------------------
%                               BULLETPROOF
%--------------------------------------------------------------------------

tol = tol^2;

[n_ref,n_dim] = size(x_ref);
id_mbr = zeros(1,n_ref);

for i_ref = 1:n_ref
    
    L2 = (x_mbr(:,1) - x_ref(i_ref,1)).^2;
    
    for i_dim = 2:n_dim
        L2 = L2 + (x_mbr(:,i_dim) - x_ref(i_ref,i_dim)).^2;
    end
    
    i_mbr = find(L2 < tol,1,'first');
    
    if ~isempty(i_mbr)
        id_mbr(i_ref) = i_mbr;
    end
    
end

id_ref = find(id_mbr);
id_mbr = id_mbr(id_ref);

end