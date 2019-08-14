

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

y0 = 0.01;
a  = 0.10;

% y0 = 0.01;
% a  = 0.30;

cCkCrd{1}  = [ -1e-3, y0 ; a, y0 ];
nCrack = 1;

%--------------------------------------------------------------------------