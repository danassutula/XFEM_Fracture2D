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
% ADD DIRECTORIES TO MATLAB SEARCH PATH
%==========================================================================

% IMPORTANT: use forward slash ('/') in paths
% in order to be platform-safe (e.g. Mac, PC)

path_main = fileparts(which('call_main.m'));
if isempty(path_main); error('path_main?');
else addpath(path_main); end

% if ~exist('path_jobsLib','var')
%     % all-jobs' source directory
%     path_jobsLib = [path_main,'/JOBS_LIBRARY'];
% end

% if ~exist('path_jobsOut','var')
%     % all-jobs' target directory
%     path_jobsOut = [path_main,'/JOBS_RESULTS'];
% end

path_main_routines = [...
    path_main,'/Routines_Assembly;',...
    path_main,'/Routines_AuxInput;',...
    path_main,'/Routines_Crack;',...
    path_main,'/Routines_Criteria;',...
    path_main,'/Routines_Element;',...
    path_main,'/Routines_EnergyMin;',...
    path_main,'/Routines_Enrichment;',...
    path_main,'/Routines_EqnSystem;',...
    path_main,'/Routines_Gauss;',...
    path_main,'/Routines_Geometry;',...
    path_main,'/Routines_Growth;',... % path_main,'/Routines_Mesh;',...
    path_main,'/Routines_Plot;',...
    path_main,'/Routines_Post;',...
    path_main,'/Routines_Shapes;',...
    path_main,'/Routines_Topology;',... % path_main,'/other_libs;',...
    path_main,'/other_tools;'...
    ];

% addpath(genpath(path_main))
% rmpath(genpath(path_main))

addpath(path_main_routines,'-begin');
