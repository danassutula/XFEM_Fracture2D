function egnod = Topo_gamAll(elnod,eledg)
%--------------------------------------------------------------------------
%
% Topo_gamAll:
%   determines boundary (free surface) topology given an element topology
%   and element's (singlular) edge topology (e.g. [1,2;2,3;3,1]), by
%   finding edges that are not shared
%
% INPUT:
%   elnod - element topology for patch, size = [n_elems,n_elemNodes]
%   eledg - element's edge topology, size = [n_elemEdges,n_edgeNodes]
%
% OUTPUT:
%   eledg - boundary topology (element edges), size = [n_bndEdges,2]
%
%--------------------------------------------------------------------------
%                               BULLETPROOF
%--------------------------------------------------------------------------

n_elm = size(elnod,1);
[n_edg,n_egn] = size(eledg); % n. edges per element, n. nodes per edge

n_egs = n_edg*n_elm;        % n. edges (total)
egnod = zeros(n_egs,n_egn); % all edges' topo.

for i = 1:n_edg
    egnod(n_elm*(i-1)+1:n_elm*i,:) = elnod(:,eledg(i,:));
end

egnod = sortrows(sort(egnod,2));
idedg = all(diff(egnod,1)==0,2);
idedg = [~1;idedg] | [idedg;~1]; % (true for shared edges)

% get free boundary edges
egnod = egnod(~idedg,:); % (random orientation)

end