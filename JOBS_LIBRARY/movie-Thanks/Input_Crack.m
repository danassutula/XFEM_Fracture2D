

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
load([job_srcdir,'/Input_Crack/',sprintf('job_crackID(%s).mat','Thanks')]);

mTpFrz = [];

% fprintf('\n')
% fprintf('Note! Some initially frozen crack tips:\n\n')
% fprintf('\t[%2i,%2i]\n',mTpFrz)
% fprintf('\n')
% 
% fprintf('pause: ')
% for i = 4:-1:1
%    fprintf('%i ',i)
%    pause(1)
% end
% fprintf('\n\n')

% Outputs:
% nCrack
% cCkCrd
% mTpFrz

%--------------------------------------------------------------------------
