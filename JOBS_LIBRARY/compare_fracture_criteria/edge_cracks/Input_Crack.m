

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
% load([job_srcdir,'/Input_Crack/', ...
%     sprintf('job_crackID(%i).mat',job_crackID)]);

% 2-Edge cracks:

dy = 0.1;
dx = 0.3;

cCkCrd{1}  = [    -0.6, -0.5*dy ; -0.5*dx, -0.5*dy ];
cCkCrd{2}  = [  0.5*dx,  0.5*dy ;     0.6,  0.5*dy ];
nCrack = 2;

% More cracks:
%
% cCkCrd{3} = [0.1,0.2;0.1,-0.2];
% nCrack = 3;
% 
% mTpFrz = [3,1;3,2];
% 
% cCkCrd{4} = [0,-0.15;0.2,-0.15];
% nCrack = 4;
% 
% mTpFrz = [mTpFrz;4,1;4,2];


%--------------------------------------------------------------------------