function egnod = Topo_sortEdg(egnod)
%
%--------------------------------------------------------------------------
%
% Topo_sortEdg:
%
%   sorts element edges of a single boundary (wheather continious or open)
%   (e.g. as obtained by Topo_gamAll) so that they follow a consequative 
%   order, e.g. [n1,n2; n2,n3; n3,...];
%
% INPUT:
%
%   egnod = Nx2, array of element edges in random order
%
% OUTPUT:
%
%   egnod = Nx2, array of element edges in consecutive order
%
%--------------------------------------------------------------------------

ermsg = 'expected a single index: length(i) > 1';

ndunq = unique(egnod(:));
nnunq = length(ndunq);

% boundary is closed if 2*nnunq == numel(bdnod); if boundary is open, find
% the first end-edge and put it in the beginning to be the starting edge

if 2*nnunq ~= numel(egnod)
    
    if 2*(nnunq-1) ~= numel(egnod)
       error('discontinuous boundary; edges could be missing') 
    end
    
    for i = 1:length(ndunq);
        q = egnod == ndunq(i);
        
        if nnz(q) == 1
            
            [i,j] = find(q); 
            tmp = egnod(1,:); % to replace rhs
               
            if j == 1
                egnod(1,1) = egnod(i,1);
                egnod(1,2) = egnod(i,2);
            else
                egnod(1,1) = egnod(i,2);
                egnod(1,2) = egnod(i,1);
            end
            
            egnod(i,1) = tmp(1);
            egnod(i,2) = tmp(2);
            
            break
            
        end
    end
end

q = egnod; q(1,:) = 0; % starting from first edge
for k = 1:size(q,1)-1;
    
    [i,j] = find(q == egnod(k,2));
    if length(i)>1;error(ermsg);end
        
    if j == 1
        egnod(k+1,1) = q(i,1);
        egnod(k+1,2) = q(i,2);
    else
        egnod(k+1,1) = q(i,2);
        egnod(k+1,2) = q(i,1);
    end
    
    q(i,:) = 0;
    
end
end