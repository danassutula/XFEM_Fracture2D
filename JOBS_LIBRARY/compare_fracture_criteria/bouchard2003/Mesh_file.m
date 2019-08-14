
%==========================================================================
% Mesh File
%==========================================================================


%--------------------------------------------------------------------------
% What will be expected from the mesh file
%--------------------------------------------------------------------------
%
%   mNdCrd - node cooridnates;
%            size = n_nodes x 2
%   mLNodS - element connectivities/topology;
%            size = n_elements x n_element_nodes
%   vElPhz - element material/phase tag;
%            size = n_elements x 1; max = n_materials
%   cBnNod - boundary nodes;
%            cell size = n_boundaries x 1; 
%            cell member size = n_boundary_nodes x 1, or
%            cell member size = n_boundary_lines x 2
%   cBnCrd - alternative to cBnNod;
%            boundary node coordinates; 
%            cell size = n_boundaries x 1;
%            cell member size = n_boundary_nodes x 2
%
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

mesh_FileRoot = [job_srcdir,'/Input_Mesh'];
mesh_FileType = 'gmsh'; % (only gmsh version-1)

mesh_FileName = []; % please, specify below
mesh_FilePath = []; % (concat.'ed mesh file)

%--------------------------------------------------------------------------
% Specify mesh file name
%--------------------------------------------------------------------------

mesh_FileName = '2holes2crack_buchard2003_coarse.msh';
% mesh_FileName = '2holes2crack_buchard2003_medium.msh';
% mesh_FileName = '2holes2crack_buchard2003_fine.msh';

%--------------------------------------------------------------------------
% Verify existance of mesh file
%--------------------------------------------------------------------------

mesh_FilePath = [mesh_FileRoot,'/',mesh_FileName];

if ~exist(mesh_FilePath,'file')
    error('ErrorUserInput:meshFileNotFound',...
        'Unable to find mesh file.\n')
end

%--------------------------------------------------------------------------
