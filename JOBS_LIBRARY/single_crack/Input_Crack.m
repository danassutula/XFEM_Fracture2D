

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

switch job_i
    case 1
        
        % 1-Edge crack (lhs)
        cCkCrd{1}  = [-0.1, 0; 0.5, 0];
        nCrack = 1;
        
    case 2
        
        % 2-Edge cracks (dy = 0)
        cCkCrd{1}  = [-0.1, 0; 0.3, 0];
        cCkCrd{2}  = [ 1.1, 0; 0.7, 0];
        nCrack = 2;
        
    case 3
        
        % 2-Edge cracks (dy = 0.1)
        cCkCrd{1}  = [-0.1,-0.05; 0.3,-0.05];
        cCkCrd{2}  = [ 1.1, 0.05; 0.7, 0.05];
        nCrack = 2;
        
    case 4
        
        % 2-Edge cracks (dy = 0.3)
        cCkCrd{1}  = [-0.1,-0.15; 0.3,-0.15];
        cCkCrd{2}  = [ 1.1, 0.15; 0.7, 0.15];
        nCrack = 2;
end

%--------------------------------------------------------------------------
