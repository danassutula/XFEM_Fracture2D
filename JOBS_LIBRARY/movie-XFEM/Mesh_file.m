
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

% mesh_FileName = 'SquarePlate_coarse.msh'
% mesh_FileName = 'SquarePlate_medium.msh'
% mesh_FileName = 'SquarePlate_fine.msh'

% mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_COARSE.msh';
% mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_MEDIUM.msh';
% mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_FINE.msh';

% mesh_FileName = 'SquarePlate_wEdgeSlit2.msh'
% mesh_FileName = 'Plate2Holes.msh';

% switch msh_job
%     case 1
%         mesh_FileName = '2holes2crack_buchard2003_mesh=C.msh';
%     case 2
%         mesh_FileName = '2holes2crack_buchard2003_mesh=M.msh';
%     case 3
%         mesh_FileName = '2holes2crack_buchard2003_mesh=F.msh';
%      otherwise
%         error('msh_job?')
% end

% switch msh_job
%     case 1
%         mesh_FileName = '2holes2crack_buchard2003_mesh=9770.msh';
%     case 2
%         mesh_FileName = '2holes2crack_buchard2003_mesh=39112.msh';
%     case 3
%         mesh_FileName = '2holes2crack_buchard2003_mesh=155258.msh';
%      otherwise
%         error('msh_job?')
% end

% switch msh_job
%     case 1
%         mesh_FileName = 'Plate_4PointShear_no1.msh';
%     case 2
%         mesh_FileName = 'Plate_4PointShear_no2.msh';
%     case 3
%         mesh_FileName = 'Plate_4PointShear_no3.msh';
%     case 4
%         mesh_FileName = 'Plate_4PointShear_no4.msh';
%      otherwise
%         error('msh_job?')
% end

% switch msh_job
%     case 1
%         mesh_FileName = 'Plate_1EdgCrk_1Inclu_BouchardBay2003_no1.msh';
%     case 2
%         mesh_FileName = 'Plate_1EdgCrk_1Inclu_BouchardBay2003_no2.msh';
%     case 3
%         mesh_FileName = 'Plate_1EdgCrk_1Inclu_BouchardBay2003_no3.msh';
%     case 4
%         mesh_FileName = 'Plate_1EdgCrk_1Inclu_BouchardBay2003_no4.msh';
%      otherwise
%         error('msh_job?')
% end

% switch msh_job
%     case 1
%         mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_COARSE.msh';
%     case 2
%         mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_MEDIUM.msh';
%     case 3
%         mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_FINE.msh';
%     case 4
%         mesh_FileName = 'SquarePlate_2holes_moesdolbow1999_FINE2.msh';
%      otherwise
%         error('msh_job?'')
% end

% switch msh_job
%     case 0
%         mesh_FileName = 'PMMA_3Holes_msh0.msh';
%     case 1
%         mesh_FileName = 'PMMA_3Holes_msh1.msh';
%     case 2
%         mesh_FileName = 'PMMA_3Holes_msh2.msh';
%     case 3
%         mesh_FileName = 'PMMA_3Holes_msh3.msh';
%     otherwise
%         error('msh_job?')
% end
       

%--------------------------------------------------------------------------
% Verify existance of mesh file
%--------------------------------------------------------------------------

mesh_FilePath = [mesh_FileRoot,'/',mesh_FileName];

if ~exist(mesh_FilePath,'file')
    error('ErrorUserInput:meshFileNotFound',...
        'Unable to find mesh file.\n')
end

%--------------------------------------------------------------------------

