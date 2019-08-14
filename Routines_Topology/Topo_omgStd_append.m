function [ndcrd,elnod] = ...
    Topo_omgStd_append(ndcrd,elnod,ndcrd_add,elnod_add,tol)

%--------------------------------------------------------------------------

% Topo_omgAdd_union
% Topo_omgAdd_ovlap

% Topo_omgAdd:
%   gets a consistent union of two input topolgies: elnod (parent) and 
%   elnod_add (child); any overlapping regions of elnod_add are discarded.

% INPUT:
%   ndcrd     - nodal coordinates   (of parent) 
%   elnod     - element topology    (of parent)
%   ndcrd_add - nodal coordiantes   (of child) 
%   elnod_add - element topology    (of child) 
%   tol       - tolerence for finding nodes

% OUTPUT:
%   ndcrd - nodal coordinates (union)
%   elnod - element topology (union)

%--------------------------------------------------------------------------

[ndmbr,ndref] = Nodes_mbr(ndcrd_add,ndcrd,tol);
ndnew = setdiff(1:size(ndcrd_add,1),ndmbr);

nnref = size(ndcrd,1); % get num. nd. (old)
ndcrd = [ndcrd;ndcrd_add(ndnew,:)]; % (new)

% get part of elnod_add that is not overlapping with elnod
elnod_add = elnod_add(~all(ismember(elnod_add,ndmbr),2),:);

for i = 1:length(ndnew) % consistent new nodes
    elnod_add(elnod_add == ndnew(i)) = nnref+i;
end

for i = 1:length(ndmbr) % o-lap. nd. re-assigned
    elnod_add(elnod_add == ndmbr(i)) = ndref(i);
end

% new topo. consistent with ndcrd
elnod = [elnod;elnod_add];

end