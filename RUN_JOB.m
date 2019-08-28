
%==========================================================================
% RUN JOB
%==========================================================================

% clear all
% close all

%--------------------------------------------------------------------------
% Include main path (to XFEM core)
%--------------------------------------------------------------------------

path_main = fileparts(which('call_main.m'));
if isempty(path_main); error('path_main?');
else addpath(path_main); end

%--------------------------------------------------------------------------
% Input job title
%--------------------------------------------------------------------------

if ~exist('job_title','var') || isempty(job_title)
   job_title = input('Enter job titile: ./JOBS_LIBRARY/','s');
end

%--------------------------------------------------------------------------
% Set core directories
%--------------------------------------------------------------------------

% IMPORTANT: use forward slash ('/') in paths
% in order to be platform-safe (e.g. UNIX, PC)

% all-jobs' source directory
path_jobsLib = [path_main,'/JOBS_LIBRARY'];

% all-jobs' target directory
path_jobsOut = [path_main,'/JOBS_RESULTS'];

% path for job-input data
job_srcdir = [path_jobsLib,'/',job_title];

% path for job-output data
job_outdir = [path_jobsOut,'/',job_title];

% n.b. a unique job sub-ID will be appended to 'job_outdir' in 'call_main'.
% n.b. for path used for saving basic output data (i.e. anything selected
% in 'job_srcdir/Input_Scope') see 'cd/Routines_AuxInput/Auxiliary_Scope'.

%--------------------------------------------------------------------------
% Get job source/target
%--------------------------------------------------------------------------

% IMPORTANT: use forward slash ('/') in paths
% in order to be platform-safe (e.g. UNIX, PC)

% executable script for running the job
job_main_m = [job_srcdir,'/JOB_MAIN.m'];

% executable script for all job's inputs
job_input_m = [job_srcdir,'/JOB_INPUT.m'];

%--------------------------------------------------------------------------
% Verify job source
%--------------------------------------------------------------------------
 
if ~exist(job_srcdir,'dir')
    error('ErrorUserInput:pathJobThisNotFound',...
    ['Unable to find current job''s directory:\n',...
    'job_srcdir = ''%s''\n'],job_srcdir)
end

if ~exist(job_main_m,'file')
    error('ErrorUserInput:scriptJobMasterNotFound',...
    ['Unable to find job''s main-script file:\n',...
    'job_main_m = ''%s''\n'],job_main_m);
end

if ~exist(job_input_m,'file')
    error('ErrorUserInput:scriptJobInputsNotFound',...
    ['Unable to find job''s input-script file:\n',...
    'job_input_m = ''%s''\n'],job_input_m);
end

%--------------------------------------------------------------------------
% Include job path
%--------------------------------------------------------------------------

% (add all paths, then remove)
addpath(genpath(path_jobsLib));
rmpath(genpath(path_jobsLib));
addpath(genpath(job_srcdir));

% % (optional)
% addpath(genpath(path_jobsOut));
% rmpath(genpath(path_jobsOut));
% addpath(genpath(job_outdir));

% Beware that generating paths using the function 'genpath()', appends a 
% semicolon (;) after each directory; it so happens to be that using a path
% that has such a termination is invalid in trying to determine whether the
% particular directory exists. For example, exist(some_existing_path,'dir') 
% == true, whereas exist(genpath(some_existing_path),'dir') == false.

%--------------------------------------------------------------------------
% Run job's master script
%--------------------------------------------------------------------------

run(job_main_m); % run([job_srcdir,'/JOB_MAIN.m'])

%--------------------------------------------------------------------------
