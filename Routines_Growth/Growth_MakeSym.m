function alpha = Growth_MakeSym(alpha,tol)
%
%--------------------------------------------------------------------------
%
% Growth_MakeSym:
%
%   does crack increment direction averaging for those increments whose
%   growth angles are almost equal in magnitude (subject to a tolerance)
%
% INPUT:
%
%   alpha = crack increment angles (any shape)
%   tol   = tolerance for comparing angle magnitudes (optional)
%
% OUTPUT:
%
%   alpha = increment angles after averaging of magnitudes
%
%--------------------------------------------------------------------------

if nargin == 1
    tol = 1e-2;
end

idnnz = find(alpha~=0);
b = abs(alpha(idnnz));
n = length(b);

isnew = true(n,1);

for i = 1:n; idsym = find(abs((b(i)-b)/b(i)) < tol & isnew);
    if length(idsym) > 1; isnew(idsym)=0; idtip=idnnz(idsym);
        alpha(idtip) = sign(alpha(idtip))*mean(b(idsym));
    end
end
end