function [grd_Fe,gd2_Fe] = GrdForceElm_crack_brn(x_inc,x_elm,rmp,N_gam,W,P)
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

[n_crs,x_crs] = GeoXElem(x_inc,x_elm);

if n_crs==2 && ~all(sign(LevelSet2Line(x_elm,x_inc))>=0)
    
    % glb. mapping
    with_MapTyp = 2;

    n_nod = size(x_elm,1);
    n_gsp = length(W);
    
    p = ones(n_nod,1);
    q = 1:n_nod;
    
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
    
    N = ShapesStd_omg(n_gsp,X_gsp);
    
    if ~isempty(rmp) % blending ramp
       N = N.*repmat(sum(N.*repmat(rmp,n_gsp,1),2),1,n_nod); 
    end
    
    % ds/dze (2D)
    % detJ = 0.5*s;
    
    % enr. jump (Hvi.)
    B_jmp = EnrJmp_Brn(x_gsp,x_inc(2,:));
    n_jmp = size(B_jmp,2); m=n_nod*n_jmp;
    
    grd_Fe = zeros(2,m);
    gd2_Fe = zeros(2,m);
    
    for i = 1:n_jmp
        
        w = B_jmp(:,i).*W*(0.5*s);
        N_jmp = sum(N.*w(:,p),1);
        
        grd_Fe(1,q) = dfx*N_jmp;
        grd_Fe(2,q) = dfy*N_jmp;
        
        gd2_Fe(1,q) = d2fx*N_jmp;
        gd2_Fe(2,q) = d2fy*N_jmp;
        
        q = q + n_nod;
        
    end
    
    grd_Fe = grd_Fe(:);
    gd2_Fe = gd2_Fe(:);

%     % optional (n_jmp=1)
%     Fe = [fx;fy]*N_jmp;
%     Fe = Fe(:);

else
    n = length(EnrJmp_Brn([0,0],[0,0]));
    grd_Fe = zeros(2*size(x_elm,1)*n,1);
    gd2_Fe = grd_Fe; % Fe = grd_Fe;
end
end