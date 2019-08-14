function ndcrd = MeshSmooth(ndcrd,elnod,eledg)
%--------------------------------------------------------------------------
%
% MeshSmooth (from Wikipedia):
%   Laplacian smoothing is an algorithm to smooth a polygonal mesh. For 
%   each vertex in a mesh, a new position is chosen based on local 
%   information (such as the position of neighbors) and the vertex is 
%   moved there. In the case that a mesh is topologically a rectangular 
%   grid (that is, each internal vertex is connected to four neighbors) 
%   then this operation produces the Laplacian of the mesh.
%
%   More formally, the smoothing operation may be described per-vertex as:
%   \bar{x}_i = 1/N \sum_{i=1}^N x_i, where N is the number of adjecent 
%   vertices to node i and \bar{x}_i is the new position for node i
%
% Based on SMOOTHMESH by Darren Engwirda (2007):
%
%--------------------------------------------------------------------------

nit = 20;
tol = 0.005;

n_ptn = size(ndcrd,1);                                                     % n. points
n_elm = size(elnod,1);                                                     % n. elements
n_edg = size(eledg,1);                                                     % n. edges per element
n_egs = n_edg*n_elm;                                                       % n. edges (total)

% Sparse connectivity
if n_edg == 3 % (T3)
    S = sparse(elnod(:,[1,1,2,2,3,3]),...
        elnod(:,[2,3,1,3,1,2]),1,n_ptn,n_ptn);               
else % n_edg == 4 (Q4)
    S = sparse(elnod(:,[1,1,2,2,3,3,4,4]),...
        elnod(:,[2,4,1,3,2,4,3,1]),1,n_ptn,n_ptn);
end
    
W = sum(S,2);                                                              % n. vertices at nodes
ndmov = find(W~=0);                                                        % nd. that should not move
egnod = zeros(n_egs,2);                                                    % initialize edge nodes

for i = 1:n_edg
    egnod(n_elm*(i-1)+1:n_elm*i,:) = elnod(:,eledg(i,:));
end

egnod = sortrows(sort(egnod,2));                                           % put shared edges next to each other
idedg = all(diff(egnod,1)==0,2);                                           % find shared edges
idedg = [~1;idedg] | [idedg;~1];                                           % true for all shared edges

egndB = egnod(~idedg,:);                                                   % boundary edges
egndI = egnod( idedg,:);                                                   % internal edges

egndI = egndI(1:2:end-1,:);                                                % unique edges
ndbnd = unique(egndB(:));                                                  % boundary nodes
ndmov = setdiff(ndmov,ndbnd);                                              % internal nodes

L0 = max(sqrt(sum((ndcrd(egndI(:,1),:)-ndcrd(egndI(:,2),:)).^2,2)),eps);   % Edge length

S = S(ndmov,:);
W = W(ndmov,[1,1]);

for itr = 1:nit
    
   ndcrd(ndmov,:) = (S*ndcrd)./W;                                          % Laplacian smoothing

   L = max(sqrt(sum((ndcrd(egndI(:,1),:)-ndcrd(egndI(:,2),:)).^2,2)),eps); % Edge length
   
   if norm((L-L0)./L,inf) < tol
      break
   end
   
   L0 = L;
   
end

if itr == nit
   warning('maximum number of iterations reached; mesh smoothing did not converge!');
end

end