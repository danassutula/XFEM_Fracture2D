
%==========================================================================
% Energy release rates (Gs_tht, Hs_tht)
%==========================================================================

if nTp2Up ~= length(vTp2Up)
    error('nTp2Up ~= length(vTp2Up)')
end

vGsTip = zeros(nTp2Up,1);
mHsTip = zeros(nTp2Up,nTp2Up);

if nTp2Up > 0

grd_Fg = sparse(nGlDof,nTp2Up);
% tips can't be solved for
p_nul  = false(nTp2Up,1);

for i_upd = 1:nTp2Up

if vTp2Up(i_upd)<=nCrack % target crack tip
    i_crk = vTp2Up(i_upd); i_tip = 1;   else
    i_crk = vTp2Up(i_upd)-nCrack; i_tip = 2;
end

% get ck. tip. el.
i = cBrStd_eTp{i_crk,i_tip};
if isempty(i); error(['Crack tip element is not available; ', ...
'i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]); end

% crack tip D-matrix
D = cDMatx{vElPhz(i(1))};

% init. rates
Gs_tht = 0;
Hs_tht = 0;

% grd_ui = -Kg\grd_Fi
grd_Fi = sparse(nGlDof,1);

%--------------------------------------------------------------------------
% Find elements
%--------------------------------------------------------------------------

if i_tip == 1
    x_inc = cCkCrd{i_crk}([2,1],:);
    p = cHvStd_sgm{i_crk}(:,2) == 1;
else
    x_inc = cCkCrd{i_crk}([end-1,end],:);
    p = cHvStd_sgm{i_crk}(:,1) == max(cHvStd_sgm{i_crk}(:,2));
end

% pure rot. elm. (not all yet)
vElShf = [cHvStd_elm{i_crk}(p);...
    cBrBln_eFl{i_crk,i_tip}];

% partial rot. elm. (too many)
vElMv0 = cHvStd_elm{i_crk}(~p);

p = ismember(mLNodS(vElStd,:),mLNodS(vElShf,:)); 
q = find(any(p,2)); p = p(q,:); % mbr. el. (std.)

% partial rot. elm.
vElMov = vElStd(q)';
tmp = mLNodS(vElMov,:);
vNdMov = unique(tmp(p));

% pure rot. elm. (all)
vElShf = unique([ vElShf; cHvBln_elm{i_crk,i_tip}; ...
    cBrStd_eFl{i_crk,i_tip};cBrStd_eTp{i_crk,i_tip}]);

q = mLNodS(vElMv0,:);
p = ismember(q,mLNodS(vElShf,:));

vElMv0 = vElMv0(any(p,2));
vNdMv0 = unique(q(p));

p = vElEnr(~ismember(vElEnr,[vElShf;vElMv0]));
if ~any(any(ismember(mLNodS(p,:),mLNodS(vElShf,:))))
    
% if any(ismember(vElMv0,[cBrBln_eSp{i_crk,i_tip};cBrStd_eSp{i_crk,i_tip}]))
%     
%     % difficult to compute due to other enrichments; most likely
%     % branch enrichment on top of Heaviside enriched kink-element;
%     % anyway, attempt to continue; this should not be too bad
%     
%     warning(['Gs-rates cannot be computed accurately: ', ...
%     'i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]);
%     
% end
    
%--------------------------------------------------------------------------
% Get energy release rates for every element
%--------------------------------------------------------------------------

dx_crd = zeros(length(vNdMov),2);
dx_crd(:,1) = x_inc(1,2)-mNdCrd(vNdMov,2);
dx_crd(:,2) = mNdCrd(vNdMov,1)-x_inc(1,1);

dx_el0 = zeros(nLNodS,2);

for i = vElMov(:)' % (el. on the periphery)
    p = mLNodS(i,:);
    
    x_elm = mNdCrd(p,:); dx_elm = dx_el0;
    for j = 1:nLNodS; k = find(vNdMov == p(j));
        if ~isempty(k); dx_elm(j,:) = dx_crd(k,:); end
    end
    
    [grd_Ke,gd2_Ke] = GrdStiffElm_std('rot', ...
     dx_elm,x_elm,D,mGsStd_omgDrv,vGsStd_omgWgt);
    
    u = mNDspS(:,p);
    u = u(:);
    
    Gs_tht = Gs_tht - 0.5*u'*grd_Ke*u;
    Hs_tht = Hs_tht - 0.5*u'*gd2_Ke*u;
    
    % p=mNDofS(:,p);
    p=2*p([1,1],:);
    p(1,:)=p(1,:)-1;
    p=p(:);
    
    grd_Fi(p) = grd_Fi(p) + grd_Ke*u;
    
end

dx_crd = zeros(length(vNdMv0),2);
dx_crd(:,1) = x_inc(1,2)-mNdCrd(vNdMv0,2);
dx_crd(:,2) = mNdCrd(vNdMv0,1)-x_inc(1,1);

% el. force enr. dofs (Hvi.)
p_frc = 2*nLNodS+1:4*nLNodS;
grd_Fe = zeros(p_frc(end),1);
gd2_Fe = grd_Fe;

for i = vElMv0(:)' % (el. at the hinge of rotation - only Hvi. enr.)
if length(cLNodE{i})==nLNodS % assert that only Heaviside enrichment
    
    p = mLNodS(i,:);
    x_elm = mNdCrd(p,:);
    dx_elm = dx_el0;
    
    for j = 1:nLNodS; k = find(vNdMv0 == p(j));
        if ~isempty(k); dx_elm(j,:) = dx_crd(k,:); end
    end
    
    if isempty(cXsElm{i,1})
        error(['No crack-element intersections had been stored; i_elm = ',...
        num2str(i),', i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]);
    end
    
    x_sub = [ cXsElm{i,1} ; x_elm ];
    x_sub = x_sub(Nodes_unq(x_sub,tol_abs),:);
    [x_sub,l_sub] = GeoSubTri_Delaun(x_sub,cXsElm{i,2});
    
    % get subcell - element member nodes
    [q,j] = Nodes_mbr(x_sub,x_elm,tol_abs);
    
    % get subcell increments
    dx_sub = zeros(size(x_sub));
    dx_sub(q,:) = dx_elm(j,:);
    
    % get ck. inc.-el. intersection
    [~,x] = GeoElm_xEdg(x_inc,x_elm);
    
    if ~isempty(x) % at most 2; use x(end,:)
        k=Nodes_mbr(x_sub,x(end,:),tol_abs);
        dx_sub(k,1) = x_inc(1,2)-x_sub(k,2);
        dx_sub(k,2) = x_sub(k,1)-x_inc(1,1);
    end
    
    if i_tip == 1
        % element specific crack segments' coordinates
        q = cHvStd_sgm{i_crk}(cHvStd_elm{i_crk} == i,:);
        x_crk = cCkCrd{i_crk}(1:q(2)+1,:); % begin at tip, end at last+1
    else % i_tip == 2
        % element specific crack segments' coordinates
        q = cHvStd_sgm{i_crk}(cHvStd_elm{i_crk} == i,:);
        x_crk = cCkCrd{i_crk}(q(1):end,:); % begin at first, end at tip
    end
    
    if strcmp(mesh_ElemType,'T3') % do analytically for this simple element
        
        [grd_Ke,gd2_Ke] = GrdStiffElm_hvi('rot',x_crk,dx_elm,x_elm,...
            dx_sub,x_sub,l_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt,D);
        
        [grd_Fe(p_frc),gd2_Fe(p_frc)] = GrdForceElm_crack_hvi_FD(i_tip,x_crk,x_elm,...
            find(any(dx_elm,2)),mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod(i_crk,:));
        
    else % too tedious to do analytically; use 4th order finite difference
        
        [grd_Ke,gd2_Ke] = GrdStiffElm_hvi_FD(i_tip,x_crk,dx_elm,x_elm,...
            dx_sub,x_sub,l_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt,D);
        
        [grd_Fe(p_frc),gd2_Fe(p_frc)] = GrdForceElm_crack_hvi_FD(i_tip,x_crk,x_elm,...
            find(any(dx_elm,2)),mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod(i_crk,:));
        
    end
    
    % enr. el. nd.
    q = cLNodE{i};
    
    u = [mNDspS(:,p),mNDspE(:,q)];
    u = u(:);
    
    Gs_tht = Gs_tht - u'*(grd_Ke*u*0.5 - grd_Fe);
    Hs_tht = Hs_tht - u'*(gd2_Ke*u*0.5 - gd2_Fe);
    
    % p = [mNDofS(:,p),mNDofE(:,q)];
    p = 2*p([1,1],:); p(1,:)=p(1,:)-1;
    p = [p,mNDofE(:,q)]; p=p(:);
    
    grd_Fi(p) = grd_Fi(p) + grd_Ke*u - grd_Fe;
    
else warning([ ...
   'Expected only Heaviside enrichment; i_elm = ',num2str(i),...
   ', i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]);

   % difficult to compute due to other enrichments; most likely branch 
   % enrichment on top of Heaviside enriched kink-element;
   
   % solutions: consider element "i" to be a pure rotation element; this is
   % okey if the crack tip segment is much larger than the trailing segment
   % inside the element; as a result the element undergoes mainly rotation
   
   vElShf(end+1) = i; % update pure rotation elements
   
end
end

for i = vElShf(:)'
    
    p = mLNodS(i,:);
    q = cLNodE{i};
    
    x_elm = mNdCrd(p,:);
    
    Ke = MtxStiffElm(x_elm,D,[cGsEnr_omgDvS{i}, ...
        cGsEnr_omgDvE{i}],cGsEnr_omgWgt{i});
    
    grd_Fe = zeros(size(Ke,1),1);
    gd2_Fe = grd_Fe;
    
    for j = 1:size(cLEnDt{i},2)
        if cLEnDt{i}(2,j) == i_tip % Brn. enr.
            
            p_frc = nLNodS+sum(cLEnDt{i}(3,1:j-1));
            p_frc = 2*p_frc+1:2*(p_frc+nLNodS*nEnJmp_brn);
            
            % k = find(cBrBln_eSp{i_crk,i_tip}==i);
            % r = cBrBln_rSp{i_crk,i_tip}(k,:); % blending ramp
            
            tmp=cBrBln_rSp{i_crk,i_tip}(cBrBln_eSp{i_crk,i_tip}==i,:);
            [grd_Fe(p_frc),gd2_Fe(p_frc)] = GrdForceElm_crack_brn(x_inc,...
                x_elm,tmp,mGsBrn_gamShp,vGsBrn_gamWgt,mCkLod(i_crk,:));
            
        elseif cLEnDt{i}(2,j) == 0 % Hvi. enr.
            
            p_frc = nLNodS+sum(cLEnDt{i}(3,1:j));
            p_frc = 2*(p_frc-nLNodS)+1:2*p_frc;
            
            [grd_Fe(p_frc),gd2_Fe(p_frc)] = GrdForceElm_crack_hvi(i_tip,x_inc,...
                x_elm,mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod(i_crk,:));
            
        end
    end
    
    u  = [mNDspS(:,p),mNDspE(:,q)];
    
    uT = [u(2,:);-u(1,:)];
    uT = uT(:); u = u(:); % u2 = -u;
    
    % Gs =  - 1/2 * (u'*Ke*u1+u1'*Ke*u)                 + u'*grd_Fe
    % Hs =  - 1/2 * (2*u1'*Ke*u1 + u'*Ke*u2+u2'*Ke*u)   + u'*gd2_Fe
    
    % Gs =  - 1/2 * (u'*Ke*u1+u1'*Ke*u)                 + u1'*Fe
    % Hs =  - 1/2 * (2*u1'*Ke*u1 + u'*Ke*u2+u2'*Ke*u)   - u' *Fe
    
    Gs_tht = Gs_tht - u'*(Ke*uT - grd_Fe);
    Hs_tht = Hs_tht + u'*(Ke*u  + gd2_Fe) - uT'*Ke*uT;
    
%     T1(  2:2*(n+1):n^2) = -1;
%     T1(n+1:2*(n+1):n^2) =  1;
%     T2(  1:n+1:n^2) = -1;
    
%     grd_Ke = Ke*T1;
%     grd_Ke = grd_Ke+grd_Ke';
%     gd2_Ke = 2*(T1'*Ke*T1-Ke);
    
    grd_Ke = zeros(size(Ke,1));
    grd_Ke(:,1:2:end) = -Ke(:,2:2:end);
    grd_Ke(:,2:2:end) =  Ke(:,1:2:end);
    grd_Ke = grd_Ke + grd_Ke';
    
    % p = [mNDofS(:,p),mNDofE(:,q)];
    p = 2*p([1,1],:); p(1,:) = p(1,:)-1;
    p = [p,mNDofE(:,q)]; p = p(:);
    
    grd_Fi(p) = grd_Fi(p) + grd_Ke*u - grd_Fe;
    
end

vGsTip(i_upd,    1) = Gs_tht;
mHsTip(i_upd,i_upd) = Hs_tht; % +Hs_inv (later)
grd_Fg(    :,i_upd) = grd_Fi;

else % some mov. el. are enr.

% too difficult to compute due to other enrichments
warning(['Gs-rates cannot be computed accurately: ', ...
'i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]);

% mark as no-update
p_nul(i_upd) = 1;

end
end

%--------------------------------------------------------------------------
% Modify 'mHsTip' for no-update tips
%--------------------------------------------------------------------------

if any(p_nul)
    if all(p_nul); mHsTip(1:nTp2Up+1:end) = 1; grd_Fg=[]; else
        p_nul = find(p_nul); p_nul = (p_nul-1)*nTp2Up+p_nul;
        mHsTip(p_nul) = sum(abs(diag(mHsTip)))/(nTp2Up-length(p_nul));
    end
end

%--------------------------------------------------------------------------
% Add the interactive part
%--------------------------------------------------------------------------

% The current mHsTip represents the self interaction; the term to be added
% is the cross interaction; the final mHsTip obtained by adding the two is
% prone to inaccuracy; this is because it so happens to be that either part
% is a large number of similar magnitude, but of opposite signs;

if ~isempty(grd_Fg)
    mHsTip = mHsTip + ...
        full( grd_Fg(vUnDof,:)'*(mGlStf(vUnDof,vUnDof)\grd_Fg(vUnDof,:)) );
end

%--------------------------------------------------------------------------
% Delete temporary variables
%--------------------------------------------------------------------------

clear tmp q p p_nul u uT grd_Ke grd_Fe gd2_Fe gd2_Ke grd_Fi grd_Fg

%--------------------------------------------------------------------------

end
