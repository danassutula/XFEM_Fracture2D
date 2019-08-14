function [mNdCrd,mElNod,vElPhz,cBnNod] = ... % ,vPtNod
    LoadMesh_gmsh(filename,elType)

% ASSUMING FORMAT OF GMSH version-1

mesh = load_gmsh(filename);

mNdCrd = mesh.POS(:,1:2);

if strcmp(elType,'T3')
    
    mElNod = mesh.TRIANGLES(1:mesh.nbTriangles,1:3);
    vElPhz = mesh.TRIANGLES(1:mesh.nbTriangles,4);
    
    % check tri. node order
    for i = 1:size(mElNod)
        
        x0 = mNdCrd(mElNod(i,1),:);
        v1 = mNdCrd(mElNod(i,2),:)-x0;
        v2 = mNdCrd(mElNod(i,3),:)-x0;
        
        if v1(1)*v2(2)-v1(2)*v2(1) < 0
           mElNod(i,[2,3]) = mElNod(i,[3,2]);
        end
        
    end
    
elseif strcmp(elType,'Q4')
    
    mElNod = mesh.QUADS(1:mesh.nbQuads,1:4);
    vElPhz = mesh.QUADS(1:mesh.nbQuads,5);
    
else
    error(['Element type: %s',elType,...
        ' is unavailable; valid element types are: T3, Q4'])
end

% material phases
vElPhz_unq = unique(vElPhz); % sorted
for i = 1:length(vElPhz_unq)
   vElPhz(vElPhz == vElPhz_unq(i)) = i; 
end

mBnNod = mesh.LINES(1:mesh.nbLines,1:2);
vBnPhz = mesh.LINES(1:mesh.nbLines,3);
vPhUnq = unique(vBnPhz); % unq. bnd.

nBnLin = length(vPhUnq);
nBnPnt = mesh.nbPoints;

nBnNod = nBnPnt+nBnLin;
cBnNod = cell(nBnNod,1);

%%% NODES FIRST, EDGES SECOND
% % sort points in the order they were defined in Gmsh
% [~,j] = sort(mesh.POINTS(1:nBnPnt,2)); k=0;
% 
% for i = 1:nBnPnt
%    cBnNod{i} = mesh.POINTS(j(i),1); 
% end
% 
% for i = nBnPnt+1:nBnNod; k=k+1;
%    cBnNod{i} = mBnNod(vBnPhz == vPhUnq(k),:);
% end

%%% EDGES FIRST, NODES SECOND
for i = 1:nBnLin
   cBnNod{i} = mBnNod(vBnPhz == vPhUnq(i),:);
end

% sort points in the order they were defined in Gmsh
[~,j] = sort(mesh.POINTS(1:nBnPnt,2)); k=0;

for i = nBnLin+1:nBnNod; k=k+1;
   cBnNod{i} = mesh.POINTS(j(k),1);
end
