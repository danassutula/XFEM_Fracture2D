
%==========================================================================
% Read Mesh Files
%==========================================================================


%--------------------------------------------------------------------------
% What is expected from the mesh file
%--------------------------------------------------------------------------

%   mNdCrd - node cooridnates;
%            size = n_nodes x 2
%   mLNodS - element connectivities/topology;
%            size = n_elements x n_element_nodes
%   vElPhz - element material/phase tag;
%            size = n_elements x 1; max = n_materials
%   cBCNod - boundary nodes;
%            cell size = n_boundaries x 1; 
%            cell member size = n_boundary_nodes x 1, or
%            cell member size = n_boundary_lines x 2
%   cBCCrd - alternative to cBCNod;
%            boundary node coordinates; 
%            cell size = n_boundaries x 1;
%            cell member size = n_boundary_nodes x 2

%--------------------------------------------------------------------------
% Get mesh data
%--------------------------------------------------------------------------

cBCNod = []; % init. boundary nodes
cBCCrd = []; % init. bnd. nd. coord.

switch mesh_FileType
    
    case 'gmsh' % version-1
        
        [mNdCrd,mLNodS,vElPhz,cBCNod] = ...
            LoadMesh_gmsh(mesh_FilePath,mesh_ElemType);
        
    % case 'matlab'
    
    % case 'other'
        
    otherwise
        error('Unknown mesh file type: "',mesh_FileType,'"')
end

%--------------------------------------------------------------------------
% OVERRIDE: Allow only a single material phase
%--------------------------------------------------------------------------
if 0
    warning('overriding material phases');
    vElPhz = ones(size(mLNodS,1),1);
end
%--------------------------------------------------------------------------
% OVERRIDE: Get boundary nodes for a rectangular domain
%--------------------------------------------------------------------------
if 0
    
    warning('overriding boundary nodes');
    cBCNod = cell(8,1);
    
    L = max(mNdCrd(:,1)) - min(mNdCrd(:,1));
    H = max(mNdCrd(:,2)) - min(mNdCrd(:,2)); tol=max(L,H)*1e-12;
    
    cBCNod{1} = find( abs(mNdCrd(:,2)-min(mNdCrd(:,2))) < tol );
    cBCNod{2} = find( abs(mNdCrd(:,1)-max(mNdCrd(:,1))) < tol );
    cBCNod{3} = find( abs(mNdCrd(:,2)-max(mNdCrd(:,2))) < tol );
    cBCNod{4} = find( abs(mNdCrd(:,1)-min(mNdCrd(:,1))) < tol );
    
    cBCNod{5} = cBCNod{1}( abs(mNdCrd(cBCNod{1},1)-min(mNdCrd(cBCNod{1},1))) < tol );
    cBCNod{6} = cBCNod{1}( abs(mNdCrd(cBCNod{1},1)-max(mNdCrd(cBCNod{1},1))) < tol );
    cBCNod{7} = cBCNod{3}( abs(mNdCrd(cBCNod{3},1)-max(mNdCrd(cBCNod{3},1))) < tol );
    cBCNod{8} = cBCNod{3}( abs(mNdCrd(cBCNod{3},1)-min(mNdCrd(cBCNod{3},1))) < tol );
    
end
%--------------------------------------------------------------------------
% OVERRIDE: Mesh smoothing
%--------------------------------------------------------------------------
if 0
    warning('smoothing mesh');
    mNdCrd = MeshSmooth(mNdCrd,mLNodS,[1,2;2,3;3,1]);
end
%--------------------------------------------------------------------------
