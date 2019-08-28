
%==========================================================================
% CALL_MAIN
%==========================================================================


%--------------------------------------------------------------------------
% Set how warnings are displayed
%--------------------------------------------------------------------------

warning('on') % all warings on
warning('on','backtrace') % backtrace warning
warning('on','verbose') % comprehensive warning (gives msgID)

%--------------------------------------------------------------------------
% Warninings to suppress
%--------------------------------------------------------------------------

% warning('off') % suppress all warnings
% warning('off','msgID') % suppress a specific warning

try warning('off','MATLAB:delaunayTriangulation:ConsConsSplitWarnId'); end
try warning('off','MATLAB:delaunayTriangulation:DupConsWarnId'); end
try warning('off','MATLAB:delaunayTriangulation:DupPtsWarnId'); end
try warning('off','MATLAB:delaunayTriangulation:DupPtsConsUpdatedWarnId'); end
try warning('off','MATLAB:delaunayTriangulation:LoopConsWarnId'); end

%--------------------------------------------------------------------------
% Verify current job ID
%--------------------------------------------------------------------------

if ~exist('job_title','var') || isempty(job_title)
    error('Undefined job.'); % define job in 'RUN_JOB'
end

if ~exist('job_subID','var') || isempty(job_subID)
    job_subID = ['job@',datestr(now,'yymmdd-HHMM')];
elseif ~ischar(job_subID)
    job_subID = num2str(job_subID,'%i');
end

% append a job sub-ID to job's output directory
job_outdir_subID = [job_outdir,'/',job_subID];

%--------------------------------------------------------------------------
% Include core paths
%--------------------------------------------------------------------------

% IMPORTANT: use forward slash ('/') in paths
% in order to be platform-safe (e.g. Mac-PC)

addpath(genpath([path_main,'/Routines_Assembly']));
addpath(genpath([path_main,'/Routines_AuxInput']));
addpath(genpath([path_main,'/Routines_EqnSystem']));
addpath(genpath([path_main,'/Routines_EnergyMin']));

addpath(genpath([path_main,'/Routines_Enrichment']));
addpath(genpath([path_main,'/Routines_Topology']));
addpath(genpath([path_main,'/Routines_Geometry']));
addpath(genpath([path_main,'/Routines_Element']));

addpath(genpath([path_main,'/Routines_Shapes']));
addpath(genpath([path_main,'/Routines_Gauss']));
addpath(genpath([path_main,'/Routines_Mesh']));

addpath(genpath([path_main,'/Routines_Criteria']));
addpath(genpath([path_main,'/Routines_Growth']));
addpath(genpath([path_main,'/Routines_Crack']));

addpath(genpath([path_main,'/Routines_Post']));
addpath(genpath([path_main,'/Routines_Plot']));

addpath(genpath([path_main,'/other_libs']));
addpath(genpath([path_main,'/other_tools']));


%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_pre = 0;
time_asm = 0;
time_sol = 0;
time_pos = 0;
time_tot = 0; tic

%--------------------------------------------------------------------------



% no load of previous (saved) progress
if ~exist('with_backup_load','var') || ...
   (exist('with_backup_load','var') && ~with_backup_load)



%==========================================================================
% INITIALISE
%==========================================================================


%--------------------------------------------------------------------------
% Global variables
%--------------------------------------------------------------------------

Declare_Global

%--------------------------------------------------------------------------
% Run job inputs
%--------------------------------------------------------------------------

run(job_input_m); % if ~exist(job_input_m), MATLAB throws error

%--------------------------------------------------------------------------
% Auxiliary data
%--------------------------------------------------------------------------

Auxiliary_Scope
Auxiliary_Material
Auxiliary_Crack
Auxiliary_Discrete
Auxiliary_Gauss

%--------------------------------------------------------------------------
% Gauss quadrature basic shapes and weights
%--------------------------------------------------------------------------

ShapesGas_std
ShapesGas_enr

%--------------------------------------------------------------------------
% Plot mesh
%--------------------------------------------------------------------------

if plot_Mesh
    figure(fig_msh); set(fig_msh,'OuterPosition',szfig_msh);
    PlotMesh; PlotCracks; title('Initial mesh and cracks');
end

%--------------------------------------------------------------------------
% Plot cracks
%--------------------------------------------------------------------------

if plot_Cracks || save_CracksEnd || save_CracksAll
    
    figure(fig_crk); set(fig_crk,'OuterPosition',szfig_crk);
    if plot_Domain; PlotDomain; end; PlotCracks;
    title('Initial crack distribution');
    
    if save_CracksEnd || save_CracksAll
        saveas(fig_crk,[path_savedImg,'plot_CracksInitial'],'fig');
        if ~plot_Cracks; close(fig_crk); end
    end
    
end

%--------------------------------------------------------------------------
% Save cracks (initial)
%--------------------------------------------------------------------------

if save_CracksEnd || save_CracksAll
    save([path_savedVar,'var_Crack_initial'],'cCkCrd','nCrack');
end

%--------------------------------------------------------------------------
% Check for existing crack intersections
%--------------------------------------------------------------------------

% This is important for X-intersections

% The X-intersection is converted into 2 T-intersections. In other words, 
% the intersecting crack is devided into two cracks that are deflected along
% the main crack (i.e. the crack that is being intersected). Also, if the
% crack branch is so short that tip enrichment does not make sense, the tip
% enrichment is removed, and the crack tip becomes inactive.

% (!) Can resolve two intersecting cracks,
% but not a chain of crack intersections

[nCrack,cCkCrd,mTpAct,mTpRdi,mCkJun] = ...
    GeoCrk_Intersect0_ifcross(cCkCrd,mTpAct,mTpRdi,1); % f_sif

%--------------------------------------------------------------------------
% Check for impending crack intersections
%--------------------------------------------------------------------------

% minimum-distance intersections:
% crack-crack and crack-boundary

% [cCkCrd,mTpAct,mTpRdi,mCkJun] = ...
%     GeoCrk_Intersect0_mindist(cCkCrd,mTpAct,mTpRdi,mCkJun,f_tip,f_sif);

[cCkCrd,mTpAct,mTpRdi,mCkJun] = ...
    GeoCrk_Intersect0_mindist_all(cCkCrd,mTpAct,mTpRdi,mCkJun,f_tip,f_sif);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% you could load a .mat file containing crack data here
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% warning('LOADING CRACKS')
% load('./cracks_file.mat')
% -> nCrack, cCkCrd, mCkJun, mTpRdi, mTpAct
  
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Set crack face tractions for all cracks            (ROUGH IMPLEMENTATION)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if wCkLod && size(mCkLod,1) == 1
    mCkLod = mCkLod(ones(nCrack,1),:);
end

%--------------------------------------------------------------------------
% Initialise variables
%--------------------------------------------------------------------------

iStep = 0;

fCont = 0;
fStop = 0;
fStop_log = struct('msg',[],'dbg',[]);

fGlLaw_inc = 0;
fGlLaw_dir = 0;

fIter_inc = 0; % energy minimization: flag for "incremental" iterations
iIter_inc = 0; % energy minimization: count iterations of growh increments

fIter_dir = 0; % energy minimization: flag for "directional" iterations
iIter_dir = 0; % energy minimization: count iterations of growh direction

mCk2Up = false(nCrack,2); % crack tips to update (after crack propagation)
jCk2Up = []; % matrix of sorted crack tips (in order of decreasing K_eqv)
vTp2Up = []; % vector of sorted crack tips (in order of decreasing K_eqv)

nTp2Up = -1; % temporary fix to allow initial solution of system; ...
% subsequently, if nTp2Up == 0 (meaning the system did not need to be 
% updated) the solution of displacements can be skipped for efficiency

mNoInc_grw = zeros(nCrack,2); % n. of crack tip extensions segments
mNoInc_upd = zeros(nCrack,2); % -//- need to update (diff. if remeshing)

cCkGrw = cell(0,1);  % track which crack tip grow
cCkLns = cell(0,1);  % each crack tip's length
vCkLen = zeros(0,1); % total crack length

Es = zeros(0,1); % strain energy
Pi = zeros(0,1); % potential energy
Gs = zeros(0,1); % energy release rate (Gs_h)

SIF_mode1 = cell(0,1); % ck. tip SIF : K1
SIF_mode2 = cell(0,1); % ck. tip SIF : K2
SIF_equiv = cell(0,1); % equiv. SIF  : K_eqv
SIF_ratio = cell(0,1); % SIF ratio   : K2/K1 (init. now, but do at end)
JIntegral = cell(0,1);

K1 = [];
K2 = [];
Ji = [];

sCrack = [];
sStres = [];
sState = [];
sRough = [];
sTimes = [];

mCkRfn = false(nCrack,2); % crack tips that are refined
mCkRfN = zeros(nCrack,2); % n. times refined at crack tips

if with_RfnXrs % || with_RfnInc
    
    % backup for remeshing
    mNdCrd_ref = mNdCrd;
    mLNodS_ref = mLNodS;
    nElems_ref = nElems;
    
    % el. around each tip 
    cElRfn = cell(nCrack,2); % (post-refin.)
    cElRf0 = cell(nCrack,2); % ( pre-refin.)
    
    % cracks that have grown
    mIsInc = false(nCrack,2);
    
    % n.b. mIsInc ensures that a specific crack has already propagated;
    % only then refinement and crack increment halfing can be carried out;
    % otherwise, halfing of the increment would undermine the initial crack
    
end

%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_pre = toc;
time_tot = time_pre;
time_ref = time_pre;

%--------------------------------------------------------------------------




%==========================================================================
% INITIAL ASSEMBLY
%==========================================================================


%--------------------------------------------------------------------------
% Discretization summary
%--------------------------------------------------------------------------


home

fprintf('\n')
fprintf('name : %s\n',job_title)
fprintf('case : %s\n',job_subID)
fprintf('\n')
fprintf('nNdStd = %9.3g\n',nNdStd);
fprintf('nElems = %9.3g\n',nElems);
fprintf('nPhase = %9.3g\n',nPhase);
fprintf('he_ref = %9.3g\n',he_ref);
fprintf('nCrack = %9.3g\n',nCrack);
fprintf('\n')


%--------------------------------------------------------------------------
% Assemble standard part
%--------------------------------------------------------------------------

Assemble_std

%--------------------------------------------------------------------------
% Assemble enriched part
%--------------------------------------------------------------------------

Assemble_enr

%--------------------------------------------------------------------------
% Essential restraints (std.)
%--------------------------------------------------------------------------

[vFxDof_all,vFxDof_nnz,vFxVal_nnz] = ApplyBC(cFxNod,mFxDim,mFxVal);
if isempty(vFxVal_nnz); vFxVal_nnz = sparse(vFxVal_nnz); end

% n.b. issparse( sparse(A(:, 1))*  full( 1) ) == 1 % (good)
% n.b. issparse( sparse(A(:,[]))*  full([]) ) == 0 % ( bad)
% n.b. issparse( sparse(A(:,[]))*sparse([]) ) == 1 % (good)

%--------------------------------------------------------------------------
% Assemble glb. (std. + enr.)
%--------------------------------------------------------------------------

if with_Update
    
    if nGDofE > 0
        mGStfS(nGlDof,nGlDof) = 0;
        vGFrcS(nGlDof,     1) = 0;
    end
    
    mGlStf = mGlStf + mGStfS;
    vGlFrc = vGlFrc + vGFrcS;
    
    if with_RdoStd
        clear mGStfS vGFrcS
    else
        % resize back std. for reuse
        mGStfS(nGDofS+1:end,:) = [];
        mGStfS(:,nGDofS+1:end) = [];
        vGFrcS(nGDofS+1:end,:) = [];
    end
    
end

% zero if prescribed (std.)
vGFrcS(vFxDof_all) = 0;

%--------------------------------------------------------------------------


else % exist('with_backup_load','var') && with_backup_load == 1


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Load previous progress
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ~exist('job_backupID','var') || isempty(job_backupID)
    job_backupID = input('job_backupID = ./','s');
end

if exist('with_backup_save','var') && with_backup_save
    load([job_outdir_subID,'/',job_backupID]);
    with_backup_save = 1;
else
    load([job_outdir_subID,'/',job_backupID]);
    with_backup_save = 0;
end

with_backup_load = 0; % (default: no load)

switch mesh_ElemType % safe way of retaining only the used element
    case 'T3'; rmpath(genpath([path_main,'/Routines_Element/Q4']))
    case 'Q4'; rmpath(genpath([path_main,'/Routines_Element/T3']))
    otherwise; error(['Unknown element type: ',mesh_ElemType])
end

iStep = iStep-1; % (a step is added later)

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


end % ~exist('with_backup_load','var') || ...


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Save current progress
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if exist('with_backup_save','var') && with_backup_save
    
    if ~exist('job_backupID','var') || isempty(job_backupID)
        job_backupID = input('job_backupID = ./','s');
    end
    
    job_backup_delt = 300; % (sec.) 
    job_backup_time = toc; % ref. time
    
else
    with_backup_save = 0;
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




%==========================================================================
% SOLUTION PROCESS
%==========================================================================


while iStep < nStep || fCont
home; iStep = iStep+1;

% backup the entire simulation data/progress
if with_backup_save && toc-job_backup_time > job_backup_delt % (sec.)
    save([job_outdir_subID,'/',job_backupID]); job_backup_time = toc;
end

fprintf('\n')
fprintf('name : %s\n',job_title)
fprintf('case : %s\n',job_subID)

fprintf('\n')
fprintf('step : %d / %d\n',iStep,nStep)
fprintf('time : %0.0f (sec.)\n\n',toc)

if nTp2Up ~= 0 % skip resolving the system since no update has occured ...
% this facility is used only in iterations of energy minimization w.r.t.
% competing crack increments. At the start nTp2Up == -1

%--------------------------------------------------------------------------
% Assemble global (if without updating)
%--------------------------------------------------------------------------

if ~with_Update

    if nGDofE > 0
        % resize to match enr.
        mGStfS(nGlDof,nGlDof) = 0;
        vGFrcS(nGlDof,     1) = 0;
    end
    
    % combine std. & enr.
    mGlStf = mGlStf + mGStfS;
    vGlFrc = vGlFrc + vGFrcS;
    
    % resize back std. for reuse
    mGStfS(nGDofS+1:end,:) = [];
    mGStfS(:,nGDofS+1:end) = [];
    vGFrcS(nGDofS+1:end,:) = [];

end

%--------------------------------------------------------------------------
% Essential restraints (std.) - overwrite since mesh might have coarsened
%--------------------------------------------------------------------------
% if with_RfnXrs || with_RfnInc
%     
%     Standard fixed degrees of freedom 'vFxDof_all' and 'vFxDof_nnz' are
%     updated during remeshing process in script 'AssembleUpd_TopoRfn'
%     
% end
%--------------------------------------------------------------------------
% Essential restraints (enr.)
%--------------------------------------------------------------------------

% % Brn. bln. dofs (not used)
% vBlDof_brn = mNDofE(:,cat(1,cNdBln_brn{:}));
% vBlDof_brn = vBlDof_brn(:);

% Hvi. bln. dofs
vBlDof_hvi = mNDofE(:,cat(1,cNdBln_hvi{:}));
vBlDof_hvi = full(vBlDof_hvi(:));

% Hvi. singular dofs
vSgDof_hvi = mNDofE(:,cat(1,cNdEnr_hvi{:})); vSgDof_hvi = vSgDof_hvi(:);
vSgDof_hvi = full(vSgDof_hvi(diag(mGlStf(vSgDof_hvi,vSgDof_hvi)) == 0));
nSgDof_hvi = length(vSgDof_hvi);

% if ~isempty(vSgDof_hvi) % print n. sing. dofs (std.)
%     warning('Hvi. singular DOFs detected: n = %d',nSgDof_hvi)
% end

%--------------------------------------------------------------------------
% Essential restraints (all)
%--------------------------------------------------------------------------

vKnDof = [vFxDof_all;vSgDof_std;vSgDof_hvi;vBlDof_hvi];
nKnDof = length(vKnDof);

p = full(mNDofE(:,vNdEnr));
vUnDof = [(1:nGDofS)';p(:)];

vUnDof(vKnDof) = [];
nUnDof = length(vUnDof);
clear vSgDof_hvi vBlDof_hvi

%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_tot = toc;
time_del = time_tot-time_ref;
time_asm = time_asm+time_del;
time_ref = time_tot;

fprintf('\nAssembly  : %3.0f\n\n',time_del)
fprintf('\n\tDOF = %d\n\n',nUnDof)

%--------------------------------------------------------------------------
% Solve system for unknown DOF
%--------------------------------------------------------------------------


nNdTot = nNdEnr + nNdStd;
mNdDsp = zeros(nDimes,nNdTot);
mNdDsp(vFxDof_nnz) = vFxVal_nnz;

mNdDsp(vUnDof) = mGlStf(vUnDof,vUnDof)\(vGlFrc(vUnDof) ...
    - mGlStf(vUnDof,vFxDof_nnz)*vFxVal_nnz);

mNDspS = mNdDsp(:,1:nNdStd);
mNDspE = zeros(nDimes,size(mNDofE,2));
mNDspE(:,vNdEnr) = mNdDsp(:,nNdStd+1:end);

% % coordiantes of displaced nodes
% mNdDfm = mNdCrd+mNDspS';


%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_tot = toc;
time_del = time_tot-time_ref;
time_sol = time_sol+time_del;
time_ref = time_tot;

fprintf('\nSolution  : %3.0f\n\n',time_del)

%--------------------------------------------------------------------------

else
    
fprintf('\nAssembly  : BYPASSED\n\n')
fprintf('\n\tDOF = %d\n\n',nUnDof)

end % nTp2Up ~= 0




%==========================================================================
% POST-PROCESSING
%==========================================================================


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~fIter_dir % determine growth directions 
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~             
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ~fIter_inc % determine which cracks to grow
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%--------------------------------------------------------------------------
% Fracture lengths (before growth at iStep)
%--------------------------------------------------------------------------

[vCkLen(iStep,1),cCkLns{iStep,1}] = ...
    GeoCrk_Length(cCkCrd,mCkJun);
% vCkLen(iStep) = sum(cCkLns{iStep})

%--------------------------------------------------------------------------
% Save crack data
%--------------------------------------------------------------------------
if save_CracksAll
    
    sCrack.cCkCrd = cCkCrd;
    sCrack.mCkJun = mCkJun;
    sCrack.mTpRdi = mTpRdi;
    sCrack.mTpAct = mTpAct;
    
    % will be saved in sState
    % sCrack.cCkLns = cCkLns;
    % sCrack.vCkLen = vCkLen;
    
    save([path_savedVar,'var_Crack_',num2str(iStep,'%05d')],'sCrack');
    
    clear sCrack
    
end
%--------------------------------------------------------------------------
% Strain energy
%--------------------------------------------------------------------------

Es(iStep,1) = 0.5*mNdDsp(:)'*(mGlStf*mNdDsp(:));
fprintf('\n\tEs  = %0.4e',Es(iStep))

Pi(iStep,1) = Es(iStep)-mNdDsp(:)'*vGlFrc;
fprintf('\n\tPi  = %0.4e\n\n',Pi(iStep))

if iStep > 1 % average energy dissipation rate per time-step
    Gs(iStep-1,1) = (Pi(iStep-1)-Pi(iStep))/(vCkLen(iStep)-vCkLen(iStep-1));
end

if isnan(Es(iStep))
    warning('Strain energy is a NaN. Simulation stopped.')
    iStep = iStep-1; break; % the current time-step fails
    % give iStep of the solution at the previous time-step
    % at iStep, you will have valid solution parameters:
    % Es(iStep), Pi, Gs, SIF_mode1, SIF_mode2, SIF_equiv 
end

%--------------------------------------------------------------------------
% Evaluate stresses at Gauss points
%--------------------------------------------------------------------------
if with_StressGPs
    
    % size(mStres_std) = [3,n_gps]; more efficient to compute this way
    [mStres_std,mGsCrd_std] = Stress_std(mNDspS,mNdCrd,vElStd,mLNodS,...
        vElPhz,nPhase,cDMatx,mPrLod,mGsStd_omgShp,mGsStd_omgDrv);
    
    % too much memory; cast to single
    mStres_std = single(mStres_std);
    mGsCrd_std = single(mGsCrd_std);
    vMises_std = Stress_vms(mStres_std);
    
    [mStres_enr,mGsCrd_enr] = Stress_enr(mNDspS,mNDspE,mNdCrd,vElEnr,mLNodS,cLNodE,...
        vElPhz,nPhase,cDMatx,mPrLod,cGsEnr_omgShS,cGsEnr_omgDvS,cGsEnr_omgDvE);
    
    mStres_enr = single(mStres_enr);
    mGsCrd_enr = single(mGsCrd_enr);
    vMises_enr = Stress_vms(mStres_enr);
    
end
%-------------------------------------------------------------------------- 
% Save stresses
%--------------------------------------------------------------------------
if save_StressAll
    
    sStres.vMises_std = vMises_std;
    sStres.vMises_enr = vMises_enr;
    sStres.mStres_std = mStres_std;
    sStres.mStres_enr = mStres_enr;
    sStres.mGsCrd_std = mGsCrd_std;
    sStres.mGsCrd_enr = mGsCrd_enr;
    
    save([path_savedVar,'var_Stress_',num2str(iStep,'%05d')],'sStres');
    
    clear sStres % (too much memory)
    
end
%--------------------------------------------------------------------------
% Evaluate displacements at Gauss points
%--------------------------------------------------------------------------
if with_DisplcGPs
    
    [mGsDsp_std,mGsCrd_std] = Disps_std(vElStd); % uses globals
    
    % too much memory; cast to single
    mGsDsp_std = single(mGsDsp_std);
    mGsCrd_std = single(mGsCrd_std);
    
    [mGsDsp_enr,mGsCrd_enr] = Disps_enr(vElEnr); % uses globals
    
    mGsDsp_enr = single(mGsDsp_enr);
    mGsCrd_enr = single(mGsCrd_enr);
    
end
%-------------------------------------------------------------------------- 
% Save displacements
%--------------------------------------------------------------------------
if save_DisplcAll
    
    sDisps.mGsDsp_std = mGsDsp_std;
    sDisps.mGsDsp_enr = mGsDsp_enr;
    sDisps.mGsCrd_std = mGsCrd_std;
    sDisps.mGsCrd_enr = mGsCrd_enr;
    
    save([path_savedVar,'var_Displc_',num2str(iStep,'%05d')],'sDisps');
    
    clear sDisps % (too much memory)
    
end
%-------------------------------------------------------------------------- 



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Compute force-displacement on boundary
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if exist('plot_ForceDisplacement','var') && plot_ForceDisplacement
    
    if iStep == 1
        plot_ForceDisplacement_BoundaryEdges = Topo_gamSub( ...
            cBCNod{plot_ForceDisplacement_IdBoundary},mLNodS,jBNodS);
        plot_ForceDisplacement_Force = zeros(nStep,2);
    end
    
    [tx,ty] = ComputeForceOnBoundary(mNdCrd,mNDspS,mLNodS,...
        plot_ForceDisplacement_BoundaryEdges,plot_ForceDisplacement_UnitNormal,...
        cDMatx,vElPhz,mGsStd_gamShp,vGsStd_gamWgt);
    plot_ForceDisplacement_Force(iStep,:) = [tx,ty];
    
end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%--------------------------------------------------------------------------
% Plot cracks
%--------------------------------------------------------------------------
if plot_Cracks
    figure(fig_crk);
    if iStep == 1 || mod(iStep,20) == 0
        set(fig_crk,'OuterPosition',szfig_crk); clf(fig_crk);
        if plot_Domain; PlotDomain; end;
    end
    PlotCracks;
    title(sprintf('Fracture process (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------
% Plot enrichment
%--------------------------------------------------------------------------
if plot_Enriched
    figure(fig_enr);
    if iStep == 1;
        set(fig_enr,'OuterPosition',szfig_enr); clf(fig_enr);
    end 
    PlotMesh; PlotEnriched; PlotCracks;
    title(sprintf('Enrichment (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------
% Plot von-Mises (power scaling)
%--------------------------------------------------------------------------
if plot_VonMises
    figure(fig_ss1);
    if iStep == 1;
        set(fig_ss1,'OuterPosition',szfig_ss1); clf(fig_ss1);
    end
    PlotVonMises_scale;
    title(sprintf('von Mises stress (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------
% Plot von-Mises (contour)
%--------------------------------------------------------------------------
if plot_VmsContr
    figure(fig_ss2);
    if iStep == 1;
        set(fig_ss2,'OuterPosition',szfig_ss2); clf(fig_ss2);
    end
    PlotVonMises_contr; 
    title(sprintf('von Mises contours (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------
% Plot deformations
%--------------------------------------------------------------------------
if plot_Deformed
    figure(fig_dfm);
    if iStep == 1;
        set(fig_dfm,'OuterPosition',szfig_dfm); clf(fig_dfm);
    end
    PlotDeformed;
    title(sprintf('Deformation (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------
% Plot displacements
%--------------------------------------------------------------------------
if plot_Displace
    figure(fig_dsp);
    if iStep == 1;
        set(fig_dsp,'OuterPosition',szfig_dsp); clf(fig_dsp);
    end
    PlotDisplace;
    title(sprintf('Displacement magnitude (iStep = %i)',iStep));
end
%--------------------------------------------------------------------------


if ~isempty(findall(0,'Type','Figure'))
    pause(0.1) % pause to properly display figures
end


%--------------------------------------------------------------------------
% STOP SIMULATION
%--------------------------------------------------------------------------

% check energy is physical
if Es(iStep) < 0
   warning('Negative strain energy. Simulation stopped.')
   iStep = iStep-1; break 
end

if iStep > 1 && abs(log10(Es(iStep)/Es(iStep-1))) > 8
    warning('Very high energy change. Simulation stopped.')
    iStep = iStep-1; break
end

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% The J-Integral
%--------------------------------------------------------------------------

% (EVAL. TIP QUANTITIES)

% set radii for SIF eval.
R_sif = mTpRdi*f_sif;
% exclude inactive tips
R_sif(~mTpAct) = 0;

if with_JIntegral
    
    Ji = IntJ(cCkCrd,mNDspS,mNDspE,...
        mPrLod,vElPhz,cDMatx,R_sif);
    
    if wCkLod
        Ji = Ji + IntJ_CrkLodHvi(cCkCrd,mNDspE,...
            nGsHvi_gam,mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod,R_sif);
        Ji = Ji + IntJ_CrkLodBrn(cCkCrd,mNDspE,...
            nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,mCkLod,R_sif);
    end
    
    if wBdLod
        warning('Body force ignored in J-integral. Can be neglected if fine mesh.')
    end
    
    JIntegral{iStep,1} = Ji;
    
end

%--------------------------------------------------------------------------
% The M-Integral
%--------------------------------------------------------------------------

K1_aux=1; K2_aux=0;
K1 = IntM(nCrack,cCkCrd,mNDspS,mNDspE,mPrLod,...
    vElPhz,cDMatx,G,kappa,R_sif,K1_aux,K2_aux);

if wCkLod
    K1 = K1 + IntM_CrkLod(nCrack,cCkCrd,...
        nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,...
        vElPhz,G,kappa,R_sif,mCkLod,K1_aux,K2_aux);
end

K1_aux=0; K2_aux=1;
K2 = IntM(nCrack,cCkCrd,mNDspS,mNDspE,mPrLod,...
    vElPhz,cDMatx,G,kappa,R_sif,K1_aux,K2_aux);

if wCkLod
    K2 = K2 + IntM_CrkLod(nCrack,cCkCrd,...
        nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,...
        vElPhz,G,kappa,R_sif,mCkLod,K1_aux,K2_aux);
end

if wBdLod
    warning('Body force ignored in SIF. Can be neglected if fine mesh.')
end

% get SIF for tip #1
for i = find(K1(:,1))'
    K1(i,1) = K1(i,1)*E_str(vElPhz(cBrStd_eTp{i,1}(1)))/2;
    K2(i,1) = K2(i,1)*E_str(vElPhz(cBrStd_eTp{i,1}(1)))/2;
end

% get SIF for tip #2
for i = find(K1(:,2))'
    K1(i,2) = K1(i,2)*E_str(vElPhz(cBrStd_eTp{i,2}(1)))/2;
    K2(i,2) = K2(i,2)*E_str(vElPhz(cBrStd_eTp{i,2}(1)))/2;
end

% (END EVAL. TIP QUANTITIES)

if strcmp(kind_GrwCrt,'symmetric')
    K1 = Growth_MakeSym(K1); % tol = 1e-2
    K2 = Growth_MakeSym(K2);
end

SIF_mode1{iStep,1} = K1;
SIF_mode2{iStep,1} = K2;
SIF_ratio{iStep,1} = K2./K1;

%--------------------------------------------------------------------------
% Crack tip growth directions
%--------------------------------------------------------------------------

% increment directions
mTpBta = zeros(nCrack,2);

% tips in opening (assuming only these tips can grow)
q_sif = K1 > 0; 

switch kind_LawDir
    case 'maxhoop' % maximum tension
        
        for i = find(q_sif(:,1))'; mTpBta(i,1) = ...
            Growth_KinkAngle_maxTension(K2(i,1)/K1(i,1));
        end
        for i = find(q_sif(:,2))'; mTpBta(i,2) = ...
            Growth_KinkAngle_maxTension(K2(i,2)/K1(i,2));
        end
        
    case 'energy' % maximum energy release rate
        
%         warning(['The explicit maximum energy dissipation criterion ', ...
%             'by Hayashi is very similar to the maximum stress criterion'])
        
        for i = find(q_sif(:,1))'; mTpBta(i,1) = ...
            Growth_KinkAngle_minEnergy_Hayashi(K2(i,1)/K1(i,1));
        end
        for i = find(q_sif(:,2))'; mTpBta(i,2) = ...
            Growth_KinkAngle_minEnergy_Hayashi(K2(i,2)/K1(i,2));
        end
        
    case 'symmetry' % local symmetry (k2=0)
     
%         warning(['The explicit local symmetry criterion ', ...
%             'by Hayashi is very similar to the maximum stress criterion'])
        
        for i = find(q_sif(:,1))'; mTpBta(i,1) = ...
            Growth_KinkAngle_localSym_Hayashi(K2(i,1)/K1(i,1));
        end
        for i = find(q_sif(:,2))'; mTpBta(i,2) = ...
            Growth_KinkAngle_localSym_Hayashi(K2(i,2)/K1(i,2));
        end
        
    otherwise
        error('Invalid crack growth direction criterion')
end

if strcmp(kind_GrwCrt,'symmetric') % average directions
    mTpBta = Growth_MakeSym(mTpBta); % tol = 1e-2;
end

%--------------------------------------------------------------------------
% Get equivalent crack tip SIF
%--------------------------------------------------------------------------

switch kind_LawCrt
    case 'tension'
        
        % maximum stress criterion (as always - traction free crack faces)
        K_eqv = (K1.*cos(0.5*mTpBta)-3*K2.*sin(0.5*mTpBta)).*cos(0.5*mTpBta).^2;
        
    case 'energy' % (good approx.)
        
        % G-criterion based on energy
        % (for traction free crack faces)
        K_eqv = sqrt(K1.^2+1.56*K2.^2);
    
    case 'J-int'
        
        % G-criterion (in-plane)
        K_eqv = sqrt(K1.^2+K2.^2);
            
    case 'eliptic'
        
        % Palaniswamy and Knauss, 1978
        % eleptic G-criterion (empirical)
        K_eqv = sqrt(K1.^2+1.5*K2.^2);
    
    case 'Hayashi'
        
        % Hayashi and Nemat-Nasser, 1981
        % SIF at an arbitrary kink angle for mixed mode loading
        [k1,k2] = Growth_GetKinkedSIF_op2(K1,K2,mTpBta);
        K_eqv = sqrt(k1.^2+k2.^2); clear k1 k2
         
    case 'Nuismer' % energy == tension
        
        % Nuismer, 1975
        % SIF at an arbitrary kink angle for mixed mode loading
        % http://link.springer.com/10.1007/BF00038891
        [k1,k2] = Growth_GetKinkedSIF_op1(K1,K2,mTpBta);
        K_eqv = sqrt(k1.^2+k2.^2); clear k1 k2
       
    otherwise
        error('Invalid crack growth onset criterion')
end

K_eqv(~q_sif) = 0; % (compression)
SIF_equiv{iStep,1} = K_eqv;

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Stop simulation if no crack growth
%--------------------------------------------------------------------------

if fCont == 1 || all(q_sif(:)==0) || ~with_Growth
    fprintf('\n\tNo crack growth. Simulation stopped.\n\n');
    iStep = iStep-1; 
    nStep = iStep;
    break
end

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Determine cracks to grow
%--------------------------------------------------------------------------

switch kind_GrwCrt
    case 'maximum'
        
        % one crack tip with max K_eqv
        [~,vTp2Up] = max(K_eqv(:));
        
    case 'symmetric'
        
        % all crack tips within 1% of max K_eqv (equiv. to 2% of max Gs)
        vTp2Up = find(K_eqv > max(K_eqv(:))*0.98);
        
    case 'critical'
        
        % all tips with K_eqv > K_crt
        vTp2Up = find(K_eqv > K_crt);
        
    case 'all'
        
        % all tips with K1 > 0
        vTp2Up = find(K_eqv > 0);
        
    case 'custom'
        
        % all tips above the lower-bound (as fraction of max K_eqv)
        vTp2Up = find(K_eqv > max(K_eqv(:))*with_GrwCrt_inf);
        
    otherwise
        error('Invalid growth criterion')
end

%--------------------------------------------------------------------------
% Organize tips in discending order of energy release rate
%--------------------------------------------------------------------------

% For use in intersection management.
%
% 'vTp2Up' is different from 'mCk2Up' in that 'vTp2Up' is sorted in terms
% of decreasing K_eqv, whereas 'mCk2Up' is just a logical representation.
% The first utility of 'vTp2Up' is in crack intersection management that
% is carried out in the order of 'vTp2Up', i.e. more critical crack tips
% are assessed for intersections first.
%
% In many cases 'vTp2Up' is redundant because of 'mCk2Up'; however, 
% sometimes 'vTp2Up' is just more convenient to use than 'mCk2Up'.
%
% During iterations (incremental/directional) 'vTp2Up = find(mCk2Up)'.

mCk2Up(:) = 0;
mCk2Up(vTp2Up) = 1;
nTp2Up = length(vTp2Up);

[~,p] = sort(K_eqv(vTp2Up),'descend');
vTp2Up = vTp2Up(p);vTp2Up = vTp2Up(:);
jCk2Up = zeros(nTp2Up,2);

p = vTp2Up <= nCrack;
jCk2Up(p,1) = vTp2Up(p);
jCk2Up(p,2) = 1;

p = vTp2Up  > nCrack;
jCk2Up(p,1) = vTp2Up(p) - nCrack;
jCk2Up(p,2) = 2;

%--------------------------------------------------------------------------
% Grow cracks
%--------------------------------------------------------------------------

for i = find(mCk2Up(:,1))'
    
    n = cCkCrd{i}(1,:)-cCkCrd{i}(2,:);
    n = sqrt(n(1)^2+n(2)^2)\n;
    
    % rotate to global direction
    n = [cos(mTpBta(i,1)),sin(mTpBta(i,1))]*[n;-n(2),n(1)];
    cCkCrd{i} = [max(f_inc*mTpRdi(i,1),with_CkIncMax)*n + ...
        cCkCrd{i}(1,:); cCkCrd{i}];
    
    if mCkJun(i,2)
        mCkJun(i,2) = mCkJun(i,2) + 1;
    end
    
end

for i = find(mCk2Up(:,2))'
    
    n = cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:);
    n = sqrt(n(1)^2+n(2)^2)\n;
    
    % rotate to global direction and increment
    n = [cos(mTpBta(i,2)),sin(mTpBta(i,2))]*[n;-n(2),n(1)];
    cCkCrd{i} = [cCkCrd{i}; cCkCrd{i}(end,:) + ...
        max(f_inc*mTpRdi(i,2),with_CkIncMax)*n];
    
end

% n. crack growth increments
mNoInc_grw = double(mCk2Up);

%--------------------------------------------------------------------------
% Check if cracks are within the simulation bounding box
%--------------------------------------------------------------------------

% check if any crack tip is beyond simulation bounds; stop if true
if exist('mCkBox','var') && ~isempty(mCkBox); Growth_Bounds; end

%--------------------------------------------------------------------------
% Save crack data describing fracture state prior to any intersections 
% (used for crack growth iterations for determining which tips grow)
%--------------------------------------------------------------------------

if with_RfnXrs % backup
    
    % mCk2Up_ref = mCk2Up;
    mTpAct_ref = mTpAct;
    mTpRdi_ref = mTpRdi;
    mCkRfN_ref = mCkRfN;
    
    if with_GLwInc
        mCkRfn_ref = mCkRfn;
        mIsInc_ref = mIsInc;
    end
    
elseif with_AdpEnr && with_GLwInc
    mTpRdi_ref = mTpRdi;
end

%--------------------------------------------------------------------------
% Check for intersections between cracks
%--------------------------------------------------------------------------

% nTp2Up > 0, very concervative (for coarse meshes and close cracks)
% nTp2Up > 1, sufficient for not so coarse meshes (or sparse cracks)

if with_RfnXrs
    
    % get intersections if cracks have crossed
    [noInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_ifcross ...
        (jCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_sif, ...
        mCkRfN<nRefine_xrs);
    
    % (the condition: mCkRfN<maxRefine can allow the type of intersections
    % that may lead to rigid body modes; this is fine providing refinement
    % can be carried out so that the correspondng cracks can be drawn in)
    
else
    
    % get intersections if cracks have crossed
    [noInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_ifcross ...
        (jCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_sif);
    
end

% upd. pure growth increments
mNoInc_grw = mNoInc_grw + noInc;


% get intersections if cracks are too close: min. dist. intersections (acute)
[noInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_mindist_acu ...
    (jCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs);

% upd. pure growth increments
mNoInc_grw = mNoInc_grw + noInc;

% get intersections if cracks are too close: min. dist. intersections (obtuse/all)
[noInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_mindist ...
    (jCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs);

% upd. pure growth increments
mNoInc_grw = mNoInc_grw + noInc;


if nTp2Up>1 && with_GLwDir
    
    % get intersections if TIPS are too close (n.b. delicate choice of f_xrs)
    [noInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_tipdist ...
        (mCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs);
    
    % upd. pure growth increments
    mNoInc_grw = mNoInc_grw + noInc;
    
end


if with_RfnXrs
    
    % tips that may have to be refined
    p = mTpRdi>0 & mCkRfN<nRefine_xrs & ~mCk2Up;
    
    % check for tip proximity; if too close, will refine later
    [~,mCk2Up_tmp,~,mTpAct_tmp,mTpRdi_tmp] = GeoCrk_Intersect_tipdist(...
        p,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs*1.50);
    
    p = find(mCk2Up_tmp & ~mCk2Up);
    
    mCk2Up(p) = mCk2Up_tmp(p);
    mTpAct(p) = mTpAct_tmp(p);
    mTpRdi(p) = mTpRdi_tmp(p);
    
end


if with_BndXrs % checks for crack to domain boundary intersections
    if with_RfnXrs && with_BndXrs_freeze % if doing any remeshing
        
        % freeze crack tips if too close to domain boundary (BC)
        [noInc,cCkCrd,mTpAct] = GeoCrk_Intersect_boundary_frz ...
            (mCk2Up,cCkCrd,mTpAct,mTpRdi.*(2.^mCkRfN),f_rfn/f_tip);
        
        for i = find(noInc(:,1)<0)'; mCk2Up(i,1)=0;
            warning('Crack #%d at tip #1 has been frozen',i)
        end
        for i = find(noInc(:,2)<0)'; mCk2Up(i,2)=0;
            warning('Crack #%d at tip #2 has been frozen',i)
        end
        
        mNoInc_grw = mNoInc_grw + noInc;
        
    end
    
    % get intersection if crack tip too close to domain boundary
    [noInc,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_boundary ...
        (mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,f_xrs); % (uses cBnNod)
    
    % upd. pure growth increments
    mNoInc_grw = mNoInc_grw + noInc;
    
    if with_RfnXrs && ~with_BndXrs_refine
        % will cause to intersect and coarsen
        mCkRfN(noInc>0) = nRefine_xrs;
    end
end

% also same sgm. to enrich
mNoInc_upd = mNoInc_grw; % (default)

% n.b. mNoInc_upd will be different if remeshing or energy min. is used;
% the above persists only if simple explicit crack growh is simulated

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Addaptive refinement of fracture advance
%--------------------------------------------------------------------------

if with_RfnInc % crack-kink dependent remeshing
    if dBetaMax_mshCors == 0 % (max refinement)
        Growth_Adaptive_max
    else
        Growth_Adaptive_inc
    end
end

if with_RfnXrs % pre-intersection remeshing
    
    % pre-intersection remeshing
    Growth_Adaptive_xrs
    
    % if intersection is inevitable, cause to remesh one last time;
    % this allows a nice coarse mesh around the crack intersection;
    % n.b. if with_RfnXrs==1, mCkRfN>0 implies mCkRfN>=nRefine_xrs
    
    mCk2Rf = mCk2Up & mTpRdi==0 & mCkRfN>0; % (mCk2Rf is local)
    if any(mCk2Rf(:)); fprintf('\n\tadaptive remeshing: intersection\n')
        mCkRfN(mCk2Rf) = 0; % (coarse remeshing since mCkRfn(mCk2Rf)==1)
    end
    
    % mark cracks that have grown
    mIsInc(logical(mNoInc_grw)) = true; 
    
end

% n.b. global mCk2Rf will be defined as: mCk2Rf = mCk2Up & mCkRfn; in the
% preceding operations (e.g. in Growth_Adaptive_inc, Growth_Adaptive_xrs)
% mCk2Rf is local and refers to a particular refinement/coarsening reason

%--------------------------------------------------------------------------



%--------------------------------------------------------------------------
% Set precedence for intersections if using an energy-based growth law
%--------------------------------------------------------------------------

% this will retain only the intersected cracks. Incremented crack tips that
% do not intersect cracks will be pulled back. The elastic field will be
% recomputed. Such that the new incremented crack configuration contains no
% crack intersections. Only then the energy minimisation can be carried out.

if with_GLwInc
    
    % get annihilated tips
    p = mCk2Up & mTpRdi == 0;
    
    if any(p(:)) % & any(mNoInc_grw(:)>0)
        % if there are such crack tips
        
        % pull back other tips
        q = mCk2Up & mTpAct & mNoInc_grw == 1;
        
        if any(q(:))
            
            mCk2Up(q) = false;
            vTp2Up = find(mCk2Up);
            nTp2Up = length(vTp2Up);
            
            for i = find(q(:,1))'; cCkCrd{i}=cCkCrd{i}(2:end,:);
                if mCkJun(i,2); mCkJun(i,2) = mCkJun(i,2)-1; end
            end
            
            for i = find(q(:,2))'
                cCkCrd{i} = cCkCrd{i}(1:end-1,:);
            end
            
            mNoInc_grw(q) = 0;
            mNoInc_upd(q) = 0;
            
            if with_RfnXrs
                mTpRdi(q) = mTpRdi_ref(q);
                mCkRfn(q) = mCkRfn_ref(q);
                mCkRfN(q) = mCkRfN_ref(q);
                mIsInc(q) = mIsInc_ref(q);
            end
        
        end
    end
end

%--------------------------------------------------------------------------
% Determine if an energy-based growth law needs to be used
%--------------------------------------------------------------------------

if with_GLwInc || with_GLwDir
    p = mCk2Up & mTpAct & mNoInc_grw == 1; % tips to iterate
    
    if any(p(:))
        
        if with_GLwInc % determine flag for using global law
            fGlLaw_inc = max(mTpBta(p).^2)>dBetaMin_iterInc^2;
            if fGlLaw_inc; fIter_inc = 1; mCk2Up_itr = p; else
                fprintf('\n\tbypassig iterations (inc.)\n'); end
        end
        
        if with_GLwDir
            if with_DirAvg; fGlLaw_dir = max(mTpBta(p).^2)>dBetaMin_iterDir^2;
                if ~fGlLaw_dir; fprintf('\n\tbypassig iterations (dir.)\n'); end
            else fGlLaw_dir = 1; end
        end
        
    end
end

%--------------------------------------------------------------------------



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else % if fIter_inc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Growth_IterateInc

%--------------------------------------------------------------------------



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % if ~fIter_inc
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
else % if fIter_dir
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%--------------------------------------------------------------------------
% Iterations: global energy minimization
%--------------------------------------------------------------------------

if plot_Cracks
    figure(fig_crk);
    PlotCracks_itr;
    title(sprintf('Fracture process (iStep = %i, iIter = %i)',iStep,iIter_dir));
    pause(0.1);
end

if strcmpi(kind_LawDir,'energy')
    % minimum energy
    Growth_IterateDir
elseif strcmpi(kind_LawDir,'maxhoop')
    % maximum hoop stress
    Growth_IterateDir_maxK1
elseif strcmpi(kind_LawDir,'symmetry')
    % minimum energy
    Growth_IterateDir_zeroK2
else
    % minimum energy
    Growth_IterateDir
end

%--------------------------------------------------------------------------



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % if ~fIter_dir
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_tot = toc;
time_del = time_tot-time_ref;
time_pos = time_pos+time_del;
time_ref = time_tot;

fprintf('\nPost-proc.: %3.0f\n\n',time_del)

%--------------------------------------------------------------------------




%==========================================================================
% UPDATE
%==========================================================================


%--------------------------------------------------------------------------
% Fracture tracking
%--------------------------------------------------------------------------

if fGlLaw_inc % (upd. continuously)
    if fIter_inc && iIter_dir == 0
        cCkGrw{iStep,1} = mCk2Up;
    end
elseif fGlLaw_dir
    if iIter_dir == 0 % (upd. before iter.)
        cCkGrw{iStep,1} = mCk2Up;
    end
else % std. law (upd./step)
    cCkGrw{iStep,1} = mCk2Up;
end

%--------------------------------------------------------------------------
% Update Kg and Fg
%--------------------------------------------------------------------------

if with_Update
    Assemble_upd
else
    Assemble_enr
end

if with_AdpEnr
    Growth_Adaptive_enr
end

if fGlLaw_inc && ~fIter_inc
   fGlLaw_inc=0; iIter_inc=0; % increment iterations stopped; reset
end

if fGlLaw_dir && ~fIter_dir
    if  iIter_dir > 0
        fGlLaw_dir=0; iIter_dir=0; % direction iterations stopped; reset
    else
        if any(mTpAct(mNoInc_grw(:)==1)) % do iterations for existing active tips
            fIter_dir=1; % iIter_dir=0;  % initialise iter; counter already reset
        end
    end
end

if fIter_dir
    iIter_dir = iIter_dir + 1;
elseif fIter_inc % && ~fIter_dir
    iIter_inc = iIter_inc + 1;
end

if fIter_dir || fIter_inc
    iStep = iStep - 1;
elseif iStep==nStep && any(mCk2Up(:))
    fCont = 1; % post-process and stop
end

if fStop == 1
   fprintf('\n\tExecution stopped because: \n\t%s\n\n',fStop_log.msg)
   save([job_outdir_subID,'/error_log'],'fStop_log'); break;
end

%--------------------------------------------------------------------------

end




%==========================================================================
% FINAL POST-PROCESSING
%==========================================================================


%--------------------------------------------------------------------------
% Save cracks (final, trimmed)
%--------------------------------------------------------------------------

% trim cracks at intersections (no bln. sgm.)
[~,~,cCkCrd_trm] = GeoCrk_Length(cCkCrd,mCkJun);

if save_CracksEnd || save_CracksAll
    save([path_savedVar,'var_Crack_final'],...
        'nCrack','cCkCrd','mCkJun','cCkCrd_trm');
end

%--------------------------------------------------------------------------
% Save state variables
%--------------------------------------------------------------------------

sState.Es = Es;
sState.Pi = Pi;
sState.Gs = Gs;

sState.SIF_mode1 = SIF_mode1;
sState.SIF_mode2 = SIF_mode2;
sState.SIF_equiv = SIF_equiv;
sState.SIF_ratio = SIF_ratio;

if with_JIntegral
    sState.JIntegral = JIntegral;
end

sState.cCkGrw = cCkGrw;
sState.vCkLen = vCkLen;
sState.cCkLns = cCkLns;

n = length(SIF_equiv); % last stable material state

sState.cracklength_pregrow = vCkLen(1:n);
sState.cracklength_posgrow = vCkLen(2:n);
sState.cracklength_growavg = 0.5*(vCkLen(1:end-1)+vCkLen(2:end));

sState.loadscale_wrt_Kcrt_pregrow = [];
sState.loadscale_wrt_Gcrt_posgrow = []; % after crack growth
sState.loadscale_wrt_Gcrt_growavg = []; % averaged over crack growth

sState.Pi_scaled_wrt_Gcrt_posgrow = [];
sState.Pi_scaled_wrt_Gcrt_growavg = [];

sState.Es_scaled_wrt_Gcrt_posgrow = [];
sState.Es_scaled_wrt_Gcrt_growavg = [];

for i = 1:n
    if ~isempty(SIF_equiv{i})
        sState.loadscale_wrt_Kcrt_pregrow(i,1) = ...
            K_crt./max(SIF_equiv{i}(:));
    end
end

if exist('G_crt','var')
    if with_JIntegral
        for i = 1:n-1
            sState.loadscale_wrt_Gcrt_posgrow(i,1) = ...
                sqrt(G_crt./max(JIntegral{i+1}(:)));
            
            sState.Pi_scaled_wrt_Gcrt_posgrow(i,1) = ...
               sState.loadscale_wrt_Gcrt_posgrow(i)^2*Pi(i+1);
            
            sState.Es_scaled_wrt_Gcrt_posgrow(i,1) = ...
               sState.loadscale_wrt_Gcrt_posgrow(i)^2*Es(i+1);
        end
    end
    if ~isempty(Gs)
        sState.loadscale_wrt_Gcrt_growavg = sqrt(G_crt./Gs);
        sState.Pi_scaled_wrt_Gcrt_growavg = (Pi(1:end-1)+Pi(2:end))*(0.5*G_crt)./Gs;
        sState.Es_scaled_wrt_Gcrt_growavg = (Es(1:end-1)+Es(2:end))*(0.5*G_crt)./Gs;
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot resultant load vs displacement
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if exist('plot_ForceDisplacement','var') && plot_ForceDisplacement
    
plot_ForceDisplacement_Force = plot_ForceDisplacement_Force(1:iStep+1,:);
sState.CritForceReaction_scaled_wrt_Gcrt_growavg = [0,0; 0.5*( ...
    plot_ForceDisplacement_Force(2:end,:)+plot_ForceDisplacement_Force(1:end-1,:))];

sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,1) = ...
    sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,1).* ...
    [1;sState.loadscale_wrt_Gcrt_growavg];

sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,2) = ...
    sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,2).* ...
    [1;sState.loadscale_wrt_Gcrt_growavg];

sState.CritLoadingFactor_scaled_wrt_Gcrt_growavg = ...
    [0;sState.loadscale_wrt_Gcrt_growavg];

figure(1001); hold on; % load vs critical displacement
plot(sState.CritLoadingFactor_scaled_wrt_Gcrt_growavg,... 
    sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,2),'-or')

ylabel('Resultant force, f_y')
xlabel('Prescribed displacement, u_y')

figure(1002); hold on; % strain energy vs critical displacement
plot([0;sState.Es_scaled_wrt_Gcrt_growavg],... 
    sState.CritForceReaction_scaled_wrt_Gcrt_growavg(:,2),'-sb')

ylabel('Strain energy, Es')
xlabel('Prescribed displacement, u_y')

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if save_StateVarb
    save([path_savedVar,'var_State'],'sState');
end

%--------------------------------------------------------------------------
% Save fracture roughness
%--------------------------------------------------------------------------

if with_Roughness
    
    [Rq,Rv,Rp,Rt,xi,yi,y0] = GeoCrk_RoughPost ...
        (cCkCrd_trm,zeros(size(mCkJun)));
    
    sRough.Rq = Rq;
    sRough.Rv = Rv;
    sRough.Rp = Rp;
    sRough.Rt = Rt;
    sRough.xi = xi;
    sRough.yi = yi;
    sRough.y0 = y0;
    
    if save_Roughness
        save([path_savedVar,'var_Roughness'],'sRough');
    end
    
    if (plot_Roughness || save_Roughness) && length(xi) > 1
        
        figure(fig_rgh); clf(fig_rgh);
        set(fig_rgh,'OuterPosition',szfig_rgh);
        PlotRough; title('Surface profile roughness')
        
        if save_Roughness
            saveas(fig_rgh,[path_savedImg,'plot_Roughness'],'fig');
            if ~plot_Roughness; close(fig_rgh); pause(1e-6); end
        end
        
    end
end

%--------------------------------------------------------------------------
% Plot cracks
%--------------------------------------------------------------------------

if plot_Cracks || save_CracksEnd || save_CracksAll
    
    figure(fig_crk); clf(fig_crk);
    set(fig_crk,'OuterPosition',szfig_crk);
    if plot_Domain; PlotDomain; end; PlotCracks;
    title('Final fracture');
    
    % save final fracture paths
    if save_CracksEnd || save_CracksAll
        saveas(fig_crk,[path_savedImg,'plot_CracksFinal'],'fig')
        if ~plot_Cracks; close(fig_crk); pause(1e-6); end
    end
    
end

%--------------------------------------------------------------------------
% Plot potential energy
%--------------------------------------------------------------------------

if plot_Potential || save_StateVarb
    
    figure(fig_PiG); clf(fig_PiG);
    set(fig_PiG,'OuterPosition',szfig_PiG);
    PlotPotential;
    
    if save_StateVarb
        saveas(fig_PiG,[path_savedImg,'plot_Potential'],'fig');
        if ~plot_Potential; close(gcf); pause(1e-6); end
    end
    
end

%--------------------------------------------------------------------------
% Plot gloabl energy dissipation rate: Gs = -dPi/da
%--------------------------------------------------------------------------

if plot_DissipGlb || save_StateVarb
    
    figure(fig_GsG); clf(fig_GsG);
    set(fig_GsG,'OuterPosition',szfig_GsG); 
    PlotPotential_dissip; 
    
    if save_StateVarb
        saveas(fig_GsG,[path_savedImg,'plot_Dissipation'],'fig');
        if ~plot_DissipGlb; close(gcf); pause(1e-6); end
    end
    
end

%--------------------------------------------------------------------------
% Plot cracks (movie)
%--------------------------------------------------------------------------

if mov_Cracks
    MovCracks;
end

%--------------------------------------------------------------------------
% Plot von Mises (movie)
%--------------------------------------------------------------------------

if mov_VonMises
    MovVonMises;
end

%--------------------------------------------------------------------------
% Plot deformation (movie)
%--------------------------------------------------------------------------

if mov_Deformed
    MovDeformation;
end

%--------------------------------------------------------------------------
% Time it
%--------------------------------------------------------------------------

time_tot = toc;
time_del = time_tot-time_ref;
time_pos = time_pos+time_del;
    
sTimes = struct( ...
    'time_pre',time_pre,...
    'time_asm',time_asm,...
    'time_pos',time_pos,...
    'time_sol',time_sol,...
    'time_tot',time_tot);
 
%--------------------------------------------------------------------------
% Done
%--------------------------------------------------------------------------

if with_Saving
    
    save([path_saved,'/job_paths'],'path_saved',...
        'path_savedImg','path_savedMov','path_savedVar',...
        'iStep','nStep');
    
    save([path_saved,'/job_logs'],'sTimes');
    
end

home

fprintf('\n')
fprintf('name : %s\n',job_title)
fprintf('case : %s\n',job_subID)
fprintf('\n')
fprintf('step : %d / %d\n',iStep,nStep)
fprintf('time : %0.0f (sec.)\n',time_tot)
fprintf('\n')
fprintf('Elapsed time summary:\n')
fprintf('\n')
fprintf('   time_pre = %3.2f\n',time_pre)
fprintf('   time_asm = %3.2f\n',time_asm)
fprintf('   time_sol = %3.2f\n',time_sol)
fprintf('   time_pos = %3.2f\n',time_pos)
fprintf('   time_tot = %3.2f\n',time_tot)
fprintf('\n')
fprintf('DONE!');
fprintf('\n')

fprintf('\n')
fprintf('RESULTS:\n');
fprintf('\n')

if exist('sState','var') && ~isempty(sState)
    fprintf('   sState\n')
end

if exist('sCrack','var') && ~isempty(sCrack)
    fprintf('   sCrack\n')
end

if exist('sStres','var') && ~isempty(sStres)
    fprintf('   sStres\n')
end

if exist('sDisps','var') && ~isempty(sDisps)
    fprintf('   sDisps\n')
end

if exist('sRough','var') && ~isempty(sRough)
    fprintf('   sRough\n')
end

if exist('sTimes','var') && ~isempty(sTimes)
    fprintf('   sTimes\n')
end

%--------------------------------------------------------------------------
