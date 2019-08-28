
%==========================================================================
% OPEN JOB SCRIPTS
%==========================================================================
%
% Description: 
%   open all job's main input scripts
%
%==========================================================================

if ~exist('job_title','var') || isempty(job_title)
   job_title = input('Enter job titile: ./JOBS_LIBRARY/','s');
end

% all-jobs' source directory
path_jobsLib = [cd,'/JOBS_LIBRARY'];

% path for job-input data
job_srcdir = [path_jobsLib,'/',job_title];

%--------------------------------------------------------------------------
% Verify job source
%--------------------------------------------------------------------------
 
if ~exist([job_srcdir,'/JOB_MAIN.m'],'file')
    error('ErrorUserInput:pathJobThisNotFound',...
    ['Unable to find current job''s directory:\n',...
    'job_srcdir = ''%s''\n'],job_srcdir)
end

%--------------------------------------------------------------------------
% Open all required job scripts (.m)
%--------------------------------------------------------------------------

open([job_srcdir,'/JOB_MAIN.m'])
open([job_srcdir,'/Input_Scope.m'])
open([job_srcdir,'/Input_Material.m'])
open([job_srcdir,'/Input_Crack.m'])
open([job_srcdir,'/Input_BC.m'])
open([job_srcdir,'/Mesh_make.m'])
open([job_srcdir,'/Mesh_read.m'])
open([job_srcdir,'/Mesh_file.m'])
open([job_srcdir,'/JOB_MAIN.m']) % show main script

% job_title = [];

%--------------------------------------------------------------------------
