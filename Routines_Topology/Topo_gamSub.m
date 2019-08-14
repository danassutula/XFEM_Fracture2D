function bnnod = Topo_gamSub(nodes,elnod,eledg)
%--------------------------------------------------------------------------
%
% Topo_gamSub:
%   determines boundary topology (not necessarily a free surface) given
%   an element topology and element's (singlular) edge topology by
%   finding edges that are not shared and whose nodes are in "nodes"
%
% INPUT:
%   nodes - boundary (free surface) nodes, size = [n_nodes,1]
%   elnod - element topology for patch, size = [n_elems,n_elemNodes]
%   eledg - element's edge topology, size = [n_elemEdges,n_edgeNodes]
%
% OUTPUT:
%   bnnod - boundary topology (element edges), size = [n_bndEdges,2]
%
%--------------------------------------------------------------------------

[n_edg,n_egn] = size(eledg); % n. edges per element, n. nodes per edge

ndmbr = ismember(elnod,nodes);
idmbr = find(sum(ndmbr,2)>=n_egn); % n_egn is the minimum req. (safer!)

n_egs = n_edg*length(idmbr); % overestimate (same el. may have 2 edges)
bnnod = zeros(n_egs,n_egn);

k = 0;

for i = idmbr(:)'
    for j = 1:n_edg
        if ndmbr(i,eledg(j,:)) % true if all are true
            k = k+1; bnnod(k,:) = elnod(i,eledg(j,:));
        end
    end
end

bnnod = bnnod(1:k,:);

end