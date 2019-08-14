% Copyright (C) 2017, Danas Sutula
%
% This file is part of the program XFEM_Fracture2D.
%
% XFEM_Fracture2D program is free software: you can redistribute it and/or 
% modify it under the terms of the GNU Lesser General Public License as 
% published by the Free Software Foundation, either version 3 of the 
% License, or (at your option) any later version.
%
% XFEM_Fracture2D is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser
% General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with XFEM_Fracture2D. If not, see <https://www.gnu.org/licenses/>.


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
