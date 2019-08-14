
%--------------------------------------------------------------------------
% Rescale
%--------------------------------------------------------------------------

H = 0.4;
L = 1;

scale_f = 0.6;
shift_x = 0;
shift_y = 0;

x_crack = cell2mat(cCkCrd(:));

x_BL = [min(x_crack(:,1)),min(x_crack(:,2))];
x_TR = [max(x_crack(:,1)),max(x_crack(:,2))];

x_center = 0.5*(x_BL+x_TR);
x_center_final = 0.5*([0,-H/2]+[L,H/2]);
scale_factor = 0.7/(x_TR(1)-x_BL(1));

for i = 1:nCrack
   
    cCkCrd{i}(:,1) = x_center_final(1) + scale_factor*(cCkCrd{i}(:,1) - x_center(1));
    cCkCrd{i}(:,2) = x_center_final(2) + scale_factor*(cCkCrd{i}(:,2) - x_center(2));
    
end

%--------------------------------------------------------------------------
% Save Cracks
%--------------------------------------------------------------------------

% save('job_crackID(XFEM)','nCrack','cCkCrd','mTpFrz')
