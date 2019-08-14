function [grd_Fe,gd2_Fe] = GrdForceElm_crack_hvi(i_tip,x_inc,x_elm,N_gam,W,P)
%--------------------------------------------------------------------------
% (!) for pure-rotation elements only
%
% d1Fe = [N] [[F]] [-l -m; -m  l] [p ; t]
% d2Fe = [N] [[F]] [ m -l; -l -m] [p ; t]
%
% [fx fy] =  [-m l; l m] [p ; t];
%
% d1Fe = [N] [[F]] [ 0 -1;  0  1] [-m l; l m] [fx ; fy]
% d2Fe = [N] [[F]] [-1  0;  0 -1] [-m l; l m] [fx ; fy]
%
%--------------------------------------------------------------------------

if i_tip == 1
    x_inc = x_inc([2,1],:);
end % else do nothing

[n_crs,x_crs] = GeoElm_xEdg(x_inc,x_elm);

if n_crs==2 && ~all(sign(LevelSet2Line(x_elm,x_inc))>=0)

    % glb. mapping
    with_MapTyp = 2;
    
    u = x_crs(2,:)-x_crs(1,:);
    s = sqrt(u(1)*u(1)+u(2)*u(2));
    u = s\u;
    
    fx = u(1)*P(2)-u(2)*P(1);
    fy = u(1)*P(1)+u(2)*P(2);
    
    dfx =-fy;
    dfy = fx; 
    
    d2fx =-fx;
    d2fy =-fy;
    
    x_gsp = N_gam*x_crs;
    
    if with_MapTyp == 1
        % local mapping of GP's
        X_gsp = N_gam*Gauss_Glb2Lcl(x_crs,x_elm);
    else
        % global mapping of GP's
        X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
    end
    
    N = ShapesStd_omg(length(W),X_gsp);
    
    % ds/dze (2D)
    % detJ = 0.5*s;
    
    % enr. jump (Hvi.)
    H_jmp = EnrJmp_Hvi(x_gsp);
    
    w = W*(H_jmp(1)*0.5*s);
    p = ones(size(x_elm,1),1);
    N_jmp = sum(N.*w(:,p),1);
    
    grd_Fe = [ dfx; dfy]*N_jmp;
    gd2_Fe = [d2fx;d2fy]*N_jmp;
    
    grd_Fe = grd_Fe(:);
    gd2_Fe = gd2_Fe(:);
    
%     % optional
%     Fe = [fx;fy]*N_jmp;
%     Fe = Fe(:);
    
else
    grd_Fe = zeros(2*size(x_elm,1),1);
    gd2_Fe = grd_Fe; % Fe = grd_Fe;
end
end