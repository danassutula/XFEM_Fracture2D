

%==========================================================================
% Crack Distribution
%==========================================================================


%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

nCrack = 0;  % number of cracks; nCrack = length(cCkCrd)
cCkCrd = []; % cell of crack coordinates: {[x1,y1;x2,y2;...];[];...}
mTpFrz = []; % crack tips to restrain (or freeze): [i_crk, i_tip;...]
mCkBox = []; % crack growth bounding box: [x_min,y_min;x_max,y_max],
             % i.e. confines cracks to smaller computational domain

%--------------------------------------------------------------------------
% Define cracks manually
%--------------------------------------------------------------------------

% % load cracks' file
% load([job_srcdir,'/Input_Crack/',...
% 	  sprintf('job_crackID(%i).mat',job_crack_loadID)]);

% Inclined crack #1
nCrack = 1;
a=0.4; b=30*pi/180; dx=0; dy=0; % dx=a*eps*1e11; dy=a*eps*1e11
cCkCrd{1} = [0.5,0e-4;0.5,0e-4] + [-a/2,0;a/2,0]*[cos(b),sin(b);-sin(b),cos(b)];
cCkCrd{1}(:,1) = cCkCrd{1}(:,1) + dx;
cCkCrd{1}(:,2) = cCkCrd{1}(:,2) + dy;

% Inclined crack #2
nCrack = nCrack+1;
a=0.4; b=-30*pi/180; dx=0; dy=0;
cCkCrd{end+1} = [0.5,0e-4;0.5,0e-4] + [-a/2,0;a/2,0]*[cos(b),sin(b);-sin(b),cos(b)];
cCkCrd{end}(:,1) = cCkCrd{end}(:,1) + dx;
cCkCrd{end}(:,2) = cCkCrd{end}(:,2) + dy;

%--------------------------------------------------------------------------
