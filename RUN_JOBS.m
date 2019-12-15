
%==========================================================================
% RUN JOBS
%==========================================================================

clear all
close all

job_titles = { ...
    % 'several_cracks/edge/vertical_tension', ...
    'several_cracks/inclusion', ...
};

for job_title_i = 1:length(job_titles)
    job_title = job_titles{job_title_i};

    RUN_JOB

end

% exit(0)
