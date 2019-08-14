function [grd_Ke,gd2_Ke] = GrdStiffElm_std(type,dx,x,D,dNdX,W)

%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================

switch type % derivative type
    case 'inc' % 2nd deriv. are zero
        d2x = zeros(size(dx));
    case 'rot' % 2nd deriv. are non zero
        d2x = [-dx(:,2),dx(:,1)];
    otherwise
        error('Unknown type of differentiation; type = ''inc'' or ''rot''')
end

nn = size(x,1);
nd = 2*nn;

grd_Ke = zeros(nd);
gd2_Ke = zeros(nd);

% Ke = zeros(nd);

B0 = zeros(3,nd);
B1 = zeros(3,nd);
B2 = zeros(3,nd);

iJ = zeros(2);

qy = 2:2:nd;
qx = qy-1;
qd = 1:2;

for w = W(:)'
    
    dNdX_i = dNdX(qd,:);
    
    J  = dNdX_i*x;
    DJ = J(1)*J(4)-J(2)*J(3);
    
    iJ(1) =  J(4)/DJ; iJ(3) = -J(3)/DJ;
    iJ(2) = -J(2)/DJ; iJ(4) =  J(1)/DJ;
    
    del_J  = dNdX_i*dx;
    dl2_J  = dNdX_i*d2x;
    
    grd_DJ = dNdX_i(1,:)'*dNdX_i(2,:);
    grd_DJ = grd_DJ-grd_DJ';
    
    del_DJ =  dx(:,1)'*grd_DJ*x(:,2) + x(:,1)'*grd_DJ*dx(:,2);  
    dl2_DJ = d2x(:,1)'*grd_DJ*x(:,2) + 2*dx(:,1)'*grd_DJ*dx(:,2) + x(:,1)'*grd_DJ*d2x(:,2);
    
    % d1_dNdx = d1(J^-1) dNdX; (assuming dX = 0)
    % d2_dNdx = d2(J^-1) dNdX;
    
    dNdx     =  iJ*dNdX_i;
    del_dNdx = -iJ*del_J*dNdx;
    dl2_dNdx = -2*iJ*del_J*del_dNdx - iJ*dl2_J*dNdx;
    
    B0(1,qx) = dNdx(1,:);
    B0(3,qy) = dNdx(1,:);
    B0(2,qy) = dNdx(2,:);
    B0(3,qx) = dNdx(2,:);
    
    B1(1,qx) = del_dNdx(1,:);
    B1(3,qy) = del_dNdx(1,:);
    B1(2,qy) = del_dNdx(2,:);
    B1(3,qx) = del_dNdx(2,:);
    
    B2(1,qx) = dl2_dNdx(1,:);
    B2(3,qy) = dl2_dNdx(1,:);
    B2(2,qy) = dl2_dNdx(2,:);
    B2(3,qx) = dl2_dNdx(2,:);
    
    grd_Ke = grd_Ke + (B1'*D*B0+B0'*D*B1)*(DJ*w) + B0'*D*B0*(del_DJ*w);
    gd2_Ke = gd2_Ke + (B2'*D*B0+2*B1'*D*B1+B0'*D*B2)*(DJ*w) + 2*(B1'*D*B0+B0'*D*B1)*(del_DJ*w) + B0'*D*B0*(dl2_DJ*w);
    
%     Ke = Ke + B0'*D*B0*(DJ*w);
    
    qd = qd+2;
    
end
end