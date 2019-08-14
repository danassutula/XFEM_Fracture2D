function [ndRef,lNods] = Topo_omgRef(ndRef,lNods)
%==========================================================================

% OUTPUT:
% ndRef(1)   = updated ref. node
% lNods(:,:) = new element topology (elements nodes are in order)

% INPUT:
% ndref(1)   = reference node on which to build the new topology
% lnods(:,:) = element topology (elements nodes are not in order)

%==========================================================================
%                               BULLETPROOF
%==========================================================================

j_unq = unique(lNods); % (already sorted)
n_unq = length(j_unq);  % get number of new nodes

for i = 1:n_unq % loop unique nodes
    lNods(lNods == j_unq(i)) = i; % replace with unique
end

lNods = lNods + ndRef; % get new element topology
ndRef = ndRef + n_unq; % update new reference node

%==========================================================================

end