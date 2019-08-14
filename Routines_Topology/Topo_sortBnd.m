function bnedg = Topo_sortBnd(egnod)
%
%--------------------------------------------------------------------------
%
% Topo_sortBnd:
%
%   takes a matrix of element edges (edges that lie on the domain boundary
%   as obtained by Topo_gamAll) and sorts them for every boundary in a 
%   consequative order e.g. [n1,n2; n2,n3; n3, ... ; nN,n1] 
%
% INPUT:
%
%   egnod = Nx2, array of element edges in random order for all boundaries
%
% OUTPUT:
%
%   bnedg = Nx1 cell of boundary edges; bnedg{i} gives element edges in 
%       consecutive order for boundary i 
%
%--------------------------------------------------------------------------

n = size(egnod,1);
egnod_srt = zeros(n,2);

% sorted boundaries
bnedg = {}; % init.

% for easy indexing
ndord = [1,2;2,1];

while 1
    
    i = find(egnod(:,1),1,'first');
    
    if isempty(i)
        break
    end
    
    k = 1;
    
    egnod_srt(k,:) = egnod(i,:); egnod(i,:) = 0;
    
    while 1
        
        [i,j] = find(egnod == egnod_srt(k,2));
        
        if isempty(i)
            break
        end
        
        k = k + 1;
        
        egnod_srt(k,:) = egnod(i,ndord(j,:)); egnod(i,:)=0;
        
    end

    bnedg{end+1,1} = egnod_srt(1:k,:);

end