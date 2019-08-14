
if ~exist('save_DisplcAll','var')
    error('PLEASE DEFINE: ''save_DisplcAll'' inside the job''s Input_Scope.')
end

%--------------------------------------------------------------------------
% Ensure compatibility of tasks
%--------------------------------------------------------------------------

if ~exist('with_StressGPs','var')
    with_StressGPs = 0;
end

if ~exist('with_DisplcGPs','var')
    with_DisplcGPs = 0;
end

if plot_VonMises || ...
   plot_VmsContr || ...
   save_StressAll

    with_StressGPs = 1;
end

if plot_Deformed || ...
   plot_Displace || ...
   save_DisplcAll

    with_DisplcGPs = 1;
end

if plot_Roughness || ...
   save_Roughness
    
    with_Roughness = 1;
end

if mov_Cracks
    save_CracksAll = 1;
end

if mov_VonMises
    with_StressGPs = 1;
    save_StressAll = 1;
end

if mov_Deformed
    with_DisplcGPs = 1;
    save_DisplcAll = 1;
end

% if saving anything
if save_CracksEnd || ...
   save_CracksAll || ...
   save_StressAll || ...
   save_DisplcAll || ...
   save_Roughness || ...
   save_StateVarb
        
    with_Saving = 1;
else
    with_Saving = 0;
end

%--------------------------------------------------------------------------
% Define directories for saving variables, plots and movies
%--------------------------------------------------------------------------

if with_Saving
    
    % job's output dir. for basic output data
    path_saved = job_outdir_subID; % ,'/basic'
    
    % make sub-dir.'s for diff. outputs
    path_savedVar = [path_saved,'/var/'];
    path_savedImg = [path_saved,'/img/'];
    path_savedMov = [path_saved,'/mov/'];
    
    % rmv. existing (old) dir.
    if exist(path_saved,'dir')
        rmdir(path_saved,'s');
    end
    
    mkdir(path_saved);
    mkdir(path_saved,'var');
    mkdir(path_saved,'img');
    mkdir(path_saved,'mov');
    
else
    
    path_saved      = [];
    path_savedVar   = [];
    path_savedImg   = [];
    path_savedMov   = [];
    
end

%--------------------------------------------------------------------------
% Define and format figures
%--------------------------------------------------------------------------

scrsz = get(0,'ScreenSize');
scrsz = scrsz([3,4])-scrsz([1,2])+1;

% left-bottom: [1,1,scrsz(1)/2,scrsz(2)/2]
% right-bottom: [scrsz(1)/2+1,1,scrsz(1)/2,scrsz(2)/2]
% right-top: [scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2]
% left-top: [1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2]

if plot_Mesh
    fig_msh = 1; szfig_msh = [1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_msh); close(fig_msh); end
end
if plot_Cracks || save_CracksEnd || save_CracksAll
    fig_crk = 2; szfig_crk = [1,1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_crk); close(fig_crk); end
end
if plot_Enriched
    fig_enr = 3; szfig_enr = [1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_enr); close(fig_enr); end
end


if plot_Deformed
    fig_dfm = 4; szfig_dfm = [scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_dfm); close(fig_dfm); end
end
if plot_Displace
    fig_dsp = 5; szfig_dsp = [scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_dsp); close(fig_dsp); end
end


if plot_VonMises
    fig_ss1 = 6; szfig_ss1 = [scrsz(1)/2+1,1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_ss1); close(fig_ss1); end
end
if plot_VmsContr
    fig_ss2 = 7; szfig_ss2 = [scrsz(1)/2+1,1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_ss2); close(fig_ss2); end
end


if plot_Potential || save_StateVarb
    fig_PiG = 8; szfig_PiG = [1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_PiG); close(fig_PiG); end
end
if plot_DissipGlb || save_StateVarb
    fig_GsG = 9; szfig_GsG = [scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_GsG); close(fig_GsG); end
end
if plot_Roughness || save_Roughness
    fig_rgh = 10; szfig_rgh = [scrsz(1)/2+1,1,scrsz(1)/2,scrsz(2)/2];
    if ishandle(fig_rgh); close(fig_rgh); end
end


if mov_Cracks
    fig_ckM = 11; szfig_ckM = [scrsz(1)/2-640,scrsz(2)/2-360,1280,720];
    if ishandle(fig_ckM); close(fig_ckM); end
end
if mov_VonMises
    fig_ssM = 12; szfig_ssM = [scrsz(1)/2-640,scrsz(2)/2-360,1280,720];
    if ishandle(fig_ssM); close(fig_ssM); end
end
if mov_Deformed
    fig_dfM = 13; szfig_dfM = [scrsz(1)/2-640,scrsz(2)/2-360,1280,720];
    if ishandle(fig_dfM); close(fig_dfM); end
end

%--------------------------------------------------------------------------
