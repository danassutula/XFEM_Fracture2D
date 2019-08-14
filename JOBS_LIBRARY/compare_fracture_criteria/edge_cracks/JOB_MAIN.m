
%==========================================================================
% JOB MAIN
%==========================================================================

close all

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

date_stamp = datestr(now,'yymmdd_HHMM');

jobs_criteria_growth     = {'tension','energy',  'energy'};
jobs_criteria_direction  = {'maxhoop','energy','symmetry'};
jobs_criteria_isImplicit = {    false,    true,      true};

jobs_plotResults_color   = {      'r',     'b',       'g'};
jobs_plotResults_marker  = {      'x',     's',       'v'};
jobs_plotResults_handle  = [];

scrsz = get(0,'ScreenSize');
scrsz = scrsz([3,4])-scrsz([1,2])+1;

%--------------------------------------------------------------------------
% Run job cases
%--------------------------------------------------------------------------

close all
    
for job_i = [1,2,3] % ,3

    %----------------------------------------------------------------------
    % Give job-case a unique ID
    %----------------------------------------------------------------------
    
    job_subID = sprintf('criteria_direction(%s)_isImplicit(%i)',...
        jobs_criteria_direction{job_i},jobs_criteria_isImplicit{job_i});
    
    %----------------------------------------------------------------------
    % Call main (XFEM) script
    %----------------------------------------------------------------------
    
    call_main
    
    %----------------------------------------------------------------------
    % Compare fracture paths by different criteria
    %----------------------------------------------------------------------
    
    if job_i == 1
        
        figure(101); hold on; axis equal;
        set(101,'OuterPosition',[1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2]);
        title('Comparison of fracture paths by different criteria')
        
        figure(102); hold on; % axis equal;
        set(102,'OuterPosition',[scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2]);
        title('Comparison of energy dissipation by different criteria')
        
    end
    
    %% 
    
    figure(101);
    PlotCracks; % -> h_plot
    
    set(h_plot,...
        'Color',jobs_plotResults_color{job_i},  ...
        'Marker',jobs_plotResults_marker{job_i}, ...
        'MarkerFaceColor','w' ...    
    )

    jobs_plotResults_handle(job_i,1) = h_plot(1);
    legend(jobs_plotResults_handle(1:job_i,1),...
        jobs_criteria_direction(1:job_i));
    
    
    figure(102);
    PlotPotential_dissip; % -> h_plot
    
    set(h_plot,...
        'Color',jobs_plotResults_color{job_i},  ...
        'Marker',jobs_plotResults_marker{job_i}, ...
        'MarkerFaceColor','w' ...    
    )

    jobs_plotResults_handle(job_i,2) = h_plot(1);
    legend(jobs_plotResults_handle(1:job_i,2),...
           jobs_criteria_direction(1:job_i));
    
   %%
   
   fprintf('\nPausing for 3 (sec.)\n\n')
   pause(3)
   
end
