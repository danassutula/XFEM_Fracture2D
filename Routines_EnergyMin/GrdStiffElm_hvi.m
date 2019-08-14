function [grd_Ke,gd2_Ke] = GrdStiffElm_hvi(type,x_crk,...
    dx_elm,x_elm,dx_sub,x_sub,l_sub,N_sub,dNdX_sub,W_sub,D)
%==========================================================================

% missing higher order derivative terms for element Q4

%==========================================================================

switch type % derivative type
    case 'inc' % 2nd deriv. are zero
        d2x_sub = zeros(size(dx_sub));
        d2x_elm = zeros(size(dx_elm));
    case 'rot' % 2nd deriv. are non zero
        d2x_sub = [-dx_sub(:,2),dx_sub(:,1)];
        d2x_elm = [-dx_elm(:,2),dx_elm(:,1)];
    otherwise
        error('Unknown type of differentiation; type = ''inc'' or ''rot''')
end

ng = length(W_sub);
nn = size(x_elm,1);
nd = 4*nn; % +Hvi.

grd_Ke = zeros(nd);
gd2_Ke = zeros(nd);

% Ke = zeros(nd);

B0 = zeros(3,nd);
B1 = zeros(3,nd);
B2 = zeros(3,nd);

iJ = zeros(2);

qy = 2:2:nd;
qx = qy-1;
qn = 1:nn;

for i = 1:size(l_sub,1)

  x_sbi =   x_sub(l_sub(i,:),:);
 dx_sbi =  dx_sub(l_sub(i,:),:);
d2x_sbi = d2x_sub(l_sub(i,:),:);

x_gsp = N_sub*x_sbi;
X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);

[H_gsp,H_nod] = EnrFun_Hvi(x_elm,x_gsp,x_crk);
[ ~,dNdX_elm] = ShapesStd_omg(ng,X_gsp);

% resize for easy mult. with shapes
H_gsp = H_gsp(:)'; H_gsp(2,:) = H_gsp;
H_gsp = H_gsp(:);

% resize to accomodate Hvi. enrichment
dNdX_elm(end,2*nn) = 0;

for j = 1:nn
    dNdX_elm(:,nn+j) = dNdX_elm(:,j).*(H_gsp-H_nod(j));
end

qd = [1,2];

for w = W_sub(:)'
    
    dNdX_sbi = dNdX_sub(qd,:);
    DJs = det(dNdX_sbi*x_sbi);
    
    grd_DJs = dNdX_sbi(1,:)'*dNdX_sbi(2,:);
    grd_DJs = grd_DJs-grd_DJs';
    
    del_DJs = dx_sbi(:,1)'*grd_DJs*x_sbi(:,2) + ...
              x_sbi(:,1)'*grd_DJs*dx_sbi(:,2);
    
    dl2_DJs = d2x_sbi(:,1)'*grd_DJs*x_sbi(:,2)   + ...
              2*dx_sbi(:,1)'*grd_DJs*dx_sbi(:,2) + ...
              x_sbi(:,1)'*grd_DJs*d2x_sbi(:,2);
    
    dNdX = dNdX_elm(qd,qn);
    
    J  = dNdX*x_elm;
    DJ = J(1)*J(4)-J(2)*J(3);
    
    iJ(1) =  J(4)/DJ; iJ(3) = -J(3)/DJ;
    iJ(2) = -J(2)/DJ; iJ(4) =  J(1)/DJ;
    
    del_J = dNdX*dx_elm;
    dl2_J = dNdX*d2x_elm;
    
    % d1_dNdx = d1(J^-1) dNdX; (assuming dX = 0)
    % d2_dNdx = d2(J^-1) dNdX;
    
    dNdx =  iJ*dNdX_elm(qd,:);
    del_dNdx = -iJ*del_J*dNdx;                          % ERROR for Q4
    dl2_dNdx = -2*iJ*del_J*del_dNdx - iJ*dl2_J*dNdx;    % ERROR for Q4
    
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
    
    grd_Ke = grd_Ke + (B1'*D*B0+B0'*D*B1)*(DJs*w) + B0'*D*B0*(del_DJs*w);
    gd2_Ke = gd2_Ke + (B2'*D*B0+2*B1'*D*B1+B0'*D*B2)*(DJs*w) + 2*(B1'*D*B0+B0'*D*B1)*(del_DJs*w) + B0'*D*B0*(dl2_DJs*w);
    
%     Ke = Ke + B0'*D*B0*(DJs*w);

    qd = qd+2;

end
end
end