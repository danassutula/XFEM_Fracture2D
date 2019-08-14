function [s_vm,s_01,s_02] = Stress_vms(s_ij)
    
%--------------------------------------------------------------------------
%   Computes Von-Misses and principal stresses
%--------------------------------------------------------------------------

t_max = sqrt((0.5*(s_ij(1,:)-s_ij(2,:))).^2 + s_ij(3,:).^2);
s_avg = 0.5*(s_ij(1,:)+s_ij(2,:));

s_01 = s_avg + t_max;
s_02 = s_avg - t_max;

s_vm = sqrt(0.5*(s_01.^2+s_02.^2+(s_01-s_02).^2))'; % size(n,1)

end