

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

% load cracks' file
load([job_srcdir,'/Input_Crack/',...
    sprintf('job_crackID(%s).mat',job_crackID)]);

% % Inclined centre crack
% a = 0.2;
% b = 10*pi/180;
% cCkCrd{1} = [0.5,0.0;0.5,0.0] + [-a/2,0;a/2,0]*[cos(b),sin(b);-sin(b),cos(b)];
% nCrack = 1;
% 
% % 1-Edge crack (lhs)
% cCkCrd{1}  = [-0.1,-1e-9; 0.1,-1e-9];
% nCrack = 1;
% 
% % 2-Edge cracks (1)
% cCkCrd{1}  = [0.40,-0.03;  -1e-4,-0.03];
% cCkCrd{2}  = [0.60, 0.03; 1+1e-4, 0.03];
% nCrack = 2;

%--------------------------------------------------------------------------
