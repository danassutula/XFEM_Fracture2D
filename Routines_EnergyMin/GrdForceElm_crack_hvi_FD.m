function [grd_Fe,gd2_Fe] = GrdForceElm_crack_hvi_FD(...
    i_tip,x_crk,x_elm,p_ndf,N_gam,W,P) % ,order
%--------------------------------------------------------------------------
% (!) for partial-rotation elements
%--------------------------------------------------------------------------

% finite difference
order = 4; % 4 or 2

% if ~exist('order','var')
%     order = 4; % default
% elseif order ~= 4 && order ~= 2
%     order = 4; % default
% end

if i_tip ~= 1 && i_tip ~= 2
    error('Bad input; i_tip = "1" or "2" ?')
end

d = 0.01*pi/180; % small rotation
T = [cos(d),sin(d);-sin(d),cos(d)];

if order == 4
    T0 = (T*T)';
else % order == 2
    T0 = T';
end

n_gsp = length(W);
n_ndf = length(p_ndf);
n_nod = size(x_elm,1);

Fe = zeros(2*n_nod,order+1);

% for easy manipulations
p_rpt = ones(n_ndf,1)*2;
q_rpt = ones(n_nod,1);

for i_dff = 1:order+1 % (order+1 == n. samples)
    
    if i_dff > 1
        x_elm(p_ndf,:) = x_crk(p_rpt,:) + ...
            (x_elm(p_ndf,:)-x_crk(p_rpt,:))*T;
        
        if i_tip == 1
            x_crk(1,:) = x_crk(2,:) + ...
                (x_crk(1,:)-x_crk(2,:))*T;
        else % i_tip == 2
            x_crk(end,:) = x_crk(end-1,:) + ...
                (x_crk(end,:)-x_crk(end-1,:))*T;
        end
    else
        x_elm(p_ndf,:) = x_crk(p_rpt,:) + ...
            (x_elm(p_ndf,:)-x_crk(p_rpt,:))*T0;
        
        if i_tip == 1
            x_crk(1,:) = x_crk(2,:) + ...
                (x_crk(1,:)-x_crk(2,:))*T0;
        else % i_tip == 2
            x_crk(end,:) = x_crk(end-1,:) + ...
                (x_crk(end,:)-x_crk(end-1,:))*T0;
        end
    end
    
    for i_crk = 1:size(x_crk,1)-1
        
        x_cki = x_crk([i_crk,i_crk+1],:); % ck. sgm.
        [n_crs,x_crs,q_egx] = GeoXElem(x_cki,x_elm);
        
        if n_crs == 2 % (segment cuts element)
            
            u = x_crs(2,:)-x_crs(1,:);
            s = sqrt(u(1)*u(1)+u(2)*u(2));
            u = s\u;
            
            f = [u(1)*P(2)-u(2)*P(1); ...
                 u(1)*P(1)+u(2)*P(2)];
            
            % GP's along crack
            x_gsp = N_gam*x_crs;
            
            % global mapping of GP's
            X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
            
            % % local mapping of GP's
            % X_gsp = N_gam*Gauss_Glb2Lcl(x_crs,x_elm);
            
            H = EnrJmp_Hvi(0); % H = const.
            N = ShapesStd_omg(n_gsp,X_gsp);
            
            Fe(:,i_dff) = Fe(:,i_dff) + reshape( ...
             f*sum(N.*W(:,q_rpt)*(H*s*0.5),1),[],1);
            
        end
    end
end

if order == 4
    grd_Fe = (Fe(:,1)+8*(Fe(:,4)-Fe(:,2))-Fe(:,5))/(12*d);
    gd2_Fe = (-Fe(:,1)+16*(Fe(:,2)+Fe(:,4))-30*Fe(:,3)-Fe(:,5))/(12*d*d);
else % order == 2
    grd_Fe = (Fe(:,3)-Fe(:,1))/(2*d);
    gd2_Fe = (Fe(:,1)-2*Fe(:,2)+Fe(:,3))/(d*d);
end
end