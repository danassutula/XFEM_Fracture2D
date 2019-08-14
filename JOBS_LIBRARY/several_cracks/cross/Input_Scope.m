
%--------------------------------------------------------------------------
% Scope of execution
%--------------------------------------------------------------------------

nStep = 10; % number of time steps

with_MeshRead = 0; % 1 - read in , 0 - generate
with_Layered  = 0; % 1 - layered domain , 0 - monolithic

if with_Layered == 1
    with_MeshRead = 0; % unexpected outcome otherwise
end

mesh_ElemType = 'Q4'; % element types: 'T3' or 'Q4'
mesh_EnriSize = 'big'; % tip enrichment radius: 'small', 'normal' or 'big'



% Crack growth criteria:
kind_LawDir = 'maxhoop'; % maxhoop; energy; symmetry;
kind_LawCrt = 'tension'; % tension; energy; J-int; eliptic; Hayashi; Nuismer;
kind_GrwCrt = 'all'; % maximum; symmetric; critical; all; custom;

if strcmpi(kind_GrwCrt,'custom')
    % lowers threshold for growth (useful if with_GLwInc == 1)
    with_GrwCrt_inf = 0.5; % (as fraction of maximum)
end



% Numerical treatment:
with_Update = 1; % efficient updating of equations
with_RdoStd = 1; % re-assemble std. stiff. mat., i.e. no recycling
with_MapTyp = 1; % (1) local mapping (faster), (2) global mapping (see note)

% n.b. global mapping should be used if cracks cut non-rectangular Q4;
% e.g.  when meshing layers of different  densities: the transition mesh
% is non-rect., hence must use with_MapTyp = 2 for elements cut by crack.

with_AdpEnr = 0; % adjust tip enrichment radius based on local mesh size



% Energy minimization:
with_GLwInc = 0; % use brute-force energy minimisation based on crack tip Gs  
with_GLwDir = 0; % use energy minimisation to determine crack growth directions
with_DirAvg = 0; % use bi-section method to get the final increment direction

if with_GLwDir
    
    nIter_dir = 5; % max iterations allowed for tip increment direction
    dBetaTol_iterDir = (0.01*pi/180) * 1; % tol. for stopping iterations
    dBetaMin_iterDir = (0.01*pi/180) * 1; % min tip kink (with_DirAvg==1)
    
end

if with_GLwInc
    dBetaMin_iterInc = (5.00*pi/180) * 0; % min tip kink for initiation
    % RmvCritr_iterInc = @(Gs) Gs < max(Gs)*0.98; % remove less favorable tips
    % RmvCritr_iterInc = @(Gs) Gs < mean(Gs)*0.98; % default criterion
end

% n.b. if with_GLwInc == 1; should use kind_GrwCrt = 'all' or, for efficiency,
% use 'custom'; this is so to have a trial set of crack increments before the
% most energetically favorable one can be determined.



% Adaptive refinement [with_Update == 1]:
with_RfnInc = 0; % addaptive re-meshing: kink dependent
with_RfnXrs = 0; % addaptive re-meshing: pre-intersection

if with_RfnInc || with_RfnXrs
    
    nRefine_inc = 1; % max times to refine crack increment/tip-mesh
    nRefine_xrs = 3; % times to refine before intersection (>=nRefine_inc)
    
    dBetaMin_mshFine = (5.0*pi/180) * 0; % refining if kink greater than
    dBetaMax_mshCors = (1.0*pi/180) * 0; % coarsening if kink smaller than
    
    % if with_RfnXrs == 0
    %     nRefine_xrs = nRefine_inc;
    % end % (will be set by default)
    
end

with_BndXrs = 1; % check for crack-boundary intersections

if with_BndXrs
    with_BndXrs_refine = 0; % keep refining untill intersection is imminent
    with_BndXrs_freeze = 1; % freeze crack tips close to domain boundary
                            % (provided refinement has already taken place)
end



% Supplementary:
with_JIntegral = 1; % crack tip energy release rate: 'Ji' 
with_Roughness = 0; % compute fracture surface roughness; ...
% physically meaningful provided fracture percolation happens



% Saving variables (in 'path_results/job_title/job_caseID/basic')
save_CracksEnd = 1; % saves initial and final crack distributions
save_CracksAll = 0; % cracks and relevant crack data at every step
save_StressAll = 0; % stress tensor at Gauss points at every step
save_DisplcAll = 0; % GP's positions and displacements at every step
save_StateVarb = 1; % state variables: Es, Pi, Gs, SIF, ...
save_Roughness = 0; % post-split fracture roughness



% Plotting (during):
plot_Mesh     = 1; % plot mesh with initial cracks
plot_Domain   = 1; % plot domain when plotting cracks
plot_Cracks   = 1; % plot cracks at every time step
plot_Enriched = 1; % enriched elements
plot_Displace = 1; % displacement contours
plot_Deformed = 0; % Gauss points after deformation
plot_VonMises = 1; % von Mises stress
plot_VmsContr = 0; % von Mises contours



% Plotting (at end):
plot_Potential = 1; % potential energy (Pi)
plot_DissipGlb = 1; % dissipation rate (-dPi/dA)
plot_Roughness = 0; % fracture surface roughness



% Movies:
mov_Cracks   = 0; % fracture evolution movie
mov_VonMises = 0; % stress evolution movie
mov_Deformed = 0; % deformation evolution movie



% Some checks
with_Debug = 0;

%--------------------------------------------------------------------------
