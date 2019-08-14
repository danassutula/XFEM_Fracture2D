function [grd_Fe,gd2_Fe] = GrdForceElm_resid_rot(dx,x,dNdX,W,f)

%--------------------------------------------------------------------------

% FE = - int B*(D*eps_0 + sig_0) dOmg = - int B*(f) dOmg

% this function is a simplification of "GrdForceElm_resid", by taking
% advantage of the fact that the element undergoes a pure rotation (no
% distortion), which can be shown to lead to: d_det(J) = 0, and a simpler
% and more efficient evaluation of d[d/dx,d/dy]

% (!!!) if the element is enriched, the enrichment must remain constant

%--------------------------------------------------------------------------

% STRESSES TRANSFORM:
% T = [  c^2,  s^2,  2*s*c; ...
%        s^2,  c^2, -2*s*c; ...
%       -s*c,  s*c,  c^2-s^2 ];

% STRAINS TRANSFORM:
% T = [  c^2,  s^2,  s*c; ...
%        s^2,  c^2, -s*c; ...
%      -2s*c, 2s*c,  c^2-s^2 ]; % <== (!)

% STRAINS TRANSFORM (rotation):
% T = [  c^2,  s^2, -s*c; ...
%        s^2,  c^2,  s*c; ...
%       2s*c,-2s*c,  c^2-s^2 ]; % <== (!)

% STRESS/STRAIN TENSOR TRANSFORM:
% s' = [c,s;-s,c]*s*[c,-s;s,c]

% rotation matrix (works for derivatives too)
% [d/dx;d/dy]' = [c, -s; s, c]*[d/dx;d/dy];

%--------------------------------------------------------------------------

f = -f; % (since RHS is negative)

d1 = [ 0,-1; 1, 0]; % d1_[d/dx;d/dy] = d1*[d/dx;d/dy];
d2 = [-1, 0; 0,-1]; % d2_[d/dx;d/dy] = d2*[d/dx;d/dy];

nJ = size(x,1);
nN = size(dNdX,2);
nD = 2*nN; % dof.

grd_Fe = zeros(nD,1);
gd2_Fe = zeros(nD,1);

B1 = zeros(2,nN);
B2 = zeros(2,nN);
iJ = zeros(2,2);

qj = 1:nJ;
qg = 1:2;

for w = W(:)'
    
    % only for the Jackobian
    dNdX_jac = dNdX(qg,qj);
    
    J  = dNdX_jac*x;
    DJ = J(1)*J(4)-J(2)*J(3);
    
    iJ(1) =  J(4)/DJ; iJ(3) = -J(3)/DJ;
    iJ(2) = -J(2)/DJ; iJ(4) =  J(1)/DJ;
    
    dNdx = iJ*dNdX(qg,:);
    % del_dNdx = d1*dNdx;
    % dl2_dNdx = d2*dNdx;
    
    B1 = B1 + d1*dNdx*(DJ*w);
    B2 = B2 + d2*dNdx*(DJ*w);
    
    qg = qg + 2;
    
end

grd_Fe(1:2:end-1) = B1(1,:)*f(1)+B1(2,:)*f(3);
grd_Fe(2:2:end)   = B1(2,:)*f(2)+B1(1,:)*f(3);

gd2_Fe(1:2:end-1) = B2(1,:)*f(1)+B2(2,:)*f(3);
gd2_Fe(2:2:end)   = B2(2,:)*f(2)+B2(1,:)*f(3);

end