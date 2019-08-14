function [grd_Fe,gd2_Fe] = GrdForceElm_resid(type,dx,x,dNdX,W,f)

%--------------------------------------------------------------------------

% FE = - int B*(D*eps_0 + sig_0) dOmg = - int B*(f) dOmg

% (!!!) if the element is enriched, the enrichment must remain constant

%--------------------------------------------------------------------------

% STRESSES TRANSFORM:
% T = [  c^2,  s^2,  2*s*c; ...
%        s^2,  c^2, -2*s*c; ...
%       -s*c,  s*c,  c^2-s^2 ];

% STRAINS TRANSFORM:
% T = [  c^2,  s^2,  s*c; ...
%        s^2,  c^2, -s*c; ...
%      -2s*c, 2s*c,  c^2-s^2 ];

% STRESS/STRAIN TENSOR TRANSFORM:
% s' = [c,s;-s,c] * s * [c,-s;s,c]

%--------------------------------------------------------------------------

f = -f; % (since RHS is negative)

switch type % derivative type
    case 'inc' % 2nd deriv. are zero
        d2x = zeros(size(dx));
    case 'rot' % 2nd deriv. are non zero
        d2x = [-dx(:,2),dx(:,1)];
    otherwise
        error('Unknown type of differentiation; type = ''inc'' or ''rot''')
end

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
    
    del_J  = dNdX_jac*dx;
    dl2_J  = dNdX_jac*d2x;
    
    grd_DJ = dNdX_jac(1,:)'*dNdX_jac(2,:);
    grd_DJ = grd_DJ-grd_DJ';
    
    del_DJ =  dx(:,1)'*grd_DJ*x(:,2) + x(:,1)'*grd_DJ*dx(:,2);
    dl2_DJ = d2x(:,1)'*grd_DJ*x(:,2) + 2*dx(:,1)'*grd_DJ*dx(:,2) + x(:,1)'*grd_DJ*d2x(:,2);
    
    % d1_dNdx = d1(J^-1) dNdX; (assuming dX = 0)
    % d2_dNdx = d2(J^-1) dNdX;
    
    dNdx     =  iJ*dNdX(qg,:);
    del_dNdx = -iJ*del_J*dNdx;
    dl2_dNdx = -2*iJ*del_J*del_dNdx - iJ*dl2_J*dNdx;
    
    B1 = B1 + del_dNdx*(DJ*w) + dNdx*(del_DJ*w);
    B2 = B2 + dl2_dNdx*(DJ*w) + 2*del_dNdx*(del_DJ*w) + dNdx*(dl2_DJ*w);
    
    qg = qg + 2;
    
end

grd_Fe(1:2:end-1) = B1(1,:)*f(1)+B1(2,:)*f(3);
grd_Fe(2:2:end)   = B1(2,:)*f(2)+B1(1,:)*f(3);

gd2_Fe(1:2:end-1) = B2(1,:)*f(1)+B2(2,:)*f(3);
gd2_Fe(2:2:end)   = B2(2,:)*f(2)+B2(1,:)*f(3);

end