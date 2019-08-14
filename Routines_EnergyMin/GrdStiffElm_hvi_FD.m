function [grd_Ke,gd2_Ke] = GrdStiffElm_hvi_FD(i_tip,x_crk,...
    dx_elm,x_elm,dx_sub,x_sub,l_sub,N_sub,dNdX_sub,W_sub,D) % ,order
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

d = 0.01*pi/180; % small rotation
T = [cos(d),sin(d);-sin(d),cos(d)];

if order == 4
    T0 = (T*T)';
else % order == 2
    T0 = T';
end

% T  is a counter-clockwise rotation matrix
% T0 is a clockwise rotation matrix (for initiating)

n_cel = size(l_sub,1); % n. sub-cells
n_gsc = length(W_sub); % n. GP per cell

n_nod = size(x_elm,1); % n. element nodes
n_gsp = n_gsc*n_cel;   % n. GP on element

x_gsp = zeros(n_gsp,2); % GP coord. on el.
W_elm = zeros(n_gsp,1); % el. GP weights

dNdX_elm = zeros(2*n_gsp,n_nod);

% % coverted to simple vector distance
% dx_sub = [dx_sub(:,2),-dx_sub(:,1)];
% dx_elm = [dx_elm(:,2),-dx_elm(:,1)];

% identify which nodes to rotate
ndrot_sub = find(any(dx_sub,2));
ndrot_elm = find(any(dx_elm,2));

% first, rotate clockwise - starting point
if i_tip == 1
    
    % pivotal point (repeated for numerical efficiency)
    pvcrd_sub = repmat(x_crk(2,:),length(ndrot_sub),1);
    pvcrd_elm = repmat(x_crk(2,:),length(ndrot_elm),1);
    
elseif i_tip == 2
    
    % pivotal point (repeated for numerical efficiency)
    pvcrd_sub = repmat(x_crk(end-1,:),length(ndrot_sub),1);
    pvcrd_elm = repmat(x_crk(end-1,:),length(ndrot_elm),1);
    
else
    error('Bad input; i_tip = "1" or "2" ?')
end

for i = 1:order+1 % (order+1 == n. samples)
    
    if i > 1 % incremental rotations 
        x_sub(ndrot_sub,:) = pvcrd_sub + ...
            (x_sub(ndrot_sub,:)-pvcrd_sub)*T;
        x_elm(ndrot_elm,:) = pvcrd_elm + ...
            (x_elm(ndrot_elm,:)-pvcrd_elm)*T;
        
        if i_tip == 1
            x_crk(1,:) = x_crk(2,:) + ...
                (x_crk(1,:)-x_crk(2,:))*T;
        else % i_tip == 2
            x_crk(end,:) = x_crk(end-1,:) + ...
                (x_crk(end,:)-x_crk(end-1,:))*T;
        end
    else % set to an initial rotation
        x_sub(ndrot_sub,:) = pvcrd_sub + ...
            (x_sub(ndrot_sub,:)-pvcrd_sub)*T0;
        x_elm(ndrot_elm,:) = pvcrd_elm + ...
            (x_elm(ndrot_elm,:)-pvcrd_elm)*T0;
        
        if i_tip == 1
            x_crk(1,:) = x_crk(2,:) + ...
                (x_crk(1,:)-x_crk(2,:))*T0;
        else % i_tip == 2
            x_crk(end,:) = x_crk(end-1,:) + ...
                (x_crk(end,:)-x_crk(end-1,:))*T0;
        end
    end
    
    q = 1:n_gsc;
    for j = 1:n_cel
        
        J = dNdX_sub*x_sub(l_sub(j,:),:);
        x_gsp(q,:) = N_sub*x_sub(l_sub(j,:),:);
        
        W_elm(q) = (J(1:2:end-1,1).*J(2:2:end,2)-...
            J(1:2:end-1,2).*J(2:2:end,1)).*W_sub; % (true area, dA)
        
        q = q + n_gsc;
        
    end
    
    X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
    [ ~,dNdX_elm] = ShapesStd_omg(n_gsp,X_gsp);
    
    J = dNdX_elm*x_elm; % el. Jackobian
    W_elm = (J(1:2:end-1,1).*J(2:2:end,2)-...
        J(1:2:end-1,2).*J(2:2:end,1)).\W_elm;
    
    [H_gsp,H_nod] = EnrFun_Hvi(x_elm,x_gsp,x_crk);
    
    % resize for easy multiplication with shapes
    H_gsp = reshape(repmat(H_gsp(:)',2,1),[],1);
    
    % resize to accomodate enrichment
    dNdX_elm(end,2*n_nod) = 0;
    
    for j = 1:n_nod
        dNdX_elm(:,n_nod+j) = dNdX_elm(:,j).*(H_gsp-H_nod(j));
    end
    
    Ke{i} = MtxStiffElm(x_elm,D,dNdX_elm,W_elm);
    
end

if order == 4
    grd_Ke = (Ke{1}+8*(Ke{4}-Ke{2})-Ke{5})/(12*d);
    gd2_Ke = (-Ke{1}+16*(Ke{2}+Ke{4})-30*Ke{3}-Ke{5})/(12*d*d);
else % order == 2
    grd_Ke = (Ke{3}-Ke{1})/(2*d);
    gd2_Ke = (Ke{1}-2*Ke{2}+Ke{3})/(d*d);
end
end