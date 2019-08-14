
%==========================================================================
% Energy release rates (Gs_inc, Hs_inc)
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
i = cBrStd_eTp{i_crk,i_tip}(1);
if isempty(i); error(['Crack tip element is not available; ', ...
'i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]); end

% ck. tip D-matrix
D = cDMatx{vElPhz(i)};

% init. rates
Gs_inc = 0;
Hs_inc = 0;

% grd_ui = -Kg\grd_Fi
grd_Fi = sparse(nGlDof,1);

%--------------------------------------------------------------------------
% Find elements
%--------------------------------------------------------------------------

if i_tip == 1
    
    x_crk = cCkCrd{i_crk}(1:2,:);
    x_inc = x_crk([2,1],:);
    
    % inc. el. (Hvi.) for end sgm.
    p = cHvStd_sgm{i_crk}(:,2) == 1;
    
else
    
    % k = max(cHvStd_sgm{i_crk}(:,2));
    k = size(cCkCrd{i_crk},1) - 1;
    
    x_crk = cCkCrd{i_crk}(k:k+1,:);
    x_inc = x_crk;
    
    % inc. el. (Hvi.) for end sgm.
    p = cHvStd_sgm{i_crk}(:,1) == k;
    
end

% inc. el. (Hvi.) for end sgm.
el_hvi = cHvStd_elm{i_crk}(p);

% [el_shf,el_mov,rm_mov] = ... % rm_mov - logical
%     Elems_tip(mNdCrd,mLNodS,x_inc(2,:),mTpRdi(i_crk,i_tip)*(1+f_inc)*0.55);
% 
% el_hvi = cHvStd_elm{i_crk}(p);
% p = ismember(el_mov,el_hvi);
% 
% el_hvi = el_mov(p);
% mn_hvi = rm_mov(p,:);
% 
% el_std = el_mov(~p);
% mn_std = rm_mov(~p,:);

el_shf = [cBrBln_eFl{i_crk,i_tip};cBrBln_eSp{i_crk,i_tip}];
mn_std = ismember(mLNodS(vElStd,:),mLNodS(el_shf,:)); 

p = any(mn_std,2);
el_std = vElStd(p);
mn_std = mn_std(p,:);

el_hvi = el_hvi(~ismember(el_hvi,[el_shf;cBrStd_eSp{i_crk,i_tip}]));
mn_hvi = ismember(mLNodS(el_hvi,:),mLNodS(el_shf,:));

p = any(mn_hvi,2);
el_hvi = el_hvi(p);
mn_hvi = mn_hvi(p,:);

if isempty(el_hvi)
    error('could not find split elements');
end

if ~any(ismember(el_std,vElEnr))

%--------------------------------------------------------------------------
% Get energy release rates
%--------------------------------------------------------------------------

u_inc = x_inc(2,:)-x_inc(1,:); % inc. dir.
u_inc = u_inc/sqrt(u_inc(1)^2+u_inc(2)^2);

dx_el0 = zeros(nLNodS,2);
for ii = 1:length(el_std)
    
    i = el_std(ii);
    p = mLNodS(i,:);
    
    x_elm = mNdCrd(p,:);
    dx_elm = dx_el0;
    
    q = mn_std(ii,:);
    dx_elm(q,1) = u_inc(1);
    dx_elm(q,2) = u_inc(2);
    
    [grd_Ke,gd2_Ke] = GrdStiffElm_std('inc', ...
        dx_elm,x_elm,D,mGsStd_omgDrv,vGsStd_omgWgt);
    
    u = mNDspS(:,p);
    u = u(:);
    
    Gs_inc = Gs_inc - 0.5*u'*grd_Ke*u;
    Hs_inc = Hs_inc - 0.5*u'*gd2_Ke*u;
    
    % p = mNDofS(:,p);
    p=2*p([1,1],:);
    p(1,:)=p(1,:)-1;
    p=p(:);
    
    grd_Fi(p) = grd_Fi(p) + grd_Ke*u;
    
end

for ii = 1:length(el_hvi)
    
    i = el_hvi(ii);
    p = mLNodS(i,:);
    
    x_elm = mNdCrd(p,:);
    dx_elm = dx_el0;
    
    q = mn_hvi(ii,:);
    dx_elm(q,1) = u_inc(1);
    dx_elm(q,2) = u_inc(2);
    
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
        k = Nodes_mbr(x_sub,x(end,:),tol_abs);
        dx_sub(k,:) = u_inc;
    end
    
    [grd_Ke,gd2_Ke] = GrdStiffElm_hvi('inc',x_crk,dx_elm,x_elm,...
        dx_sub,x_sub,l_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt,D);
    
    % enr. el. nd.
    q = cLNodE{i};
    
    if length(q) ~= nLNodS; warning([ ... % slightly less accurate
       'Expected only Heaviside enrichment; i_elm = ',num2str(i),...
       ', i_crk = ',num2str(i_crk),', i_tip = ',num2str(i_tip)]);
       
       % enr. nodes for Hvi. only
       % get position into Hvi. enr.
       k = find(cLEnDt{i}(2,:)==0);
       
       if length(k) > 1; error([... % n_Hvi > 1 (risky)
           'No more then one Heaviside enrichment; ',...
           'instead, n_Hvi = ',num2str(length(k))]);
       end
       
       % get position into Hvi. nd./dof.
       k  = sum(cLEnDt{i}(3,1:k),2);
       q = cLNodE{i}(k-nLNodS+1:k);
       
    end
    
    u = [mNDspS(:,p),mNDspE(:,q)];
    u = u(:);
    
    Gs_inc = Gs_inc - 0.5*u'*grd_Ke*u;
    Hs_inc = Hs_inc - 0.5*u'*gd2_Ke*u;
   
    % p = [mNDofS(:,p),mNDofE(:,q)];
    p = 2*p([1,1],:); p(1,:)=p(1,:)-1;
    p = [p,mNDofE(:,q)]; p = p(:);
    
    grd_Fi(p) = grd_Fi(p) + grd_Ke*u;
    
end

vGsTip(i_upd,    1) = Gs_inc;
mHsTip(i_upd,i_upd) = Hs_inc; % +Hs_inv (later)
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
    mHsTip = mHsTip + full( grd_Fg(vUnDof,:)'* ...
        (mGlStf(vUnDof,vUnDof)\grd_Fg(vUnDof,:)) );
end

%--------------------------------------------------------------------------
% Delete temporary variables
%--------------------------------------------------------------------------

% clear tmp q p p_nul u grd_Ke gd2_Ke grd_Fi grd_Fg

%--------------------------------------------------------------------------

end

warning('!!!')

if 0
%%

grd_u = mGlStf(vUnDof,vUnDof)\grd_Fg(vUnDof);

mNdDsp(vUnDof) = mNdDsp(vUnDof) + grd_u * 0.0025;

% mNdDsp = zeros(2,nNdTot);
% mNdDsp(vUnDof) = grd_u * 0.01;

mNDspS = mNdDsp(:,1:nNdStd);
mNDspE = zeros(nDimes,size(mNDofE,2));
mNDspE(:,vNdEnr) = mNdDsp(:,nNdStd+1:end);

%%
end

% figure; hold on
% PlotMesh_dsp

% Gs = KI^2/E;
% dGsda = 2K1/E*dK1da
% K1 = sqrt(Gs*E)
% dK1 = 1/2 1/sqrt(Gs*E)*E * dGsda
% dGsda = sum(Hs,2) 

%%

KI = sqrt(vGsTip*E)
dGsda = sum(mHsTip,2)
dKIda = E/2/sqrt(vGsTip(1)*E)*dGsda

vGsTip
mHsTip

figure; axis equal
PlotEnriched
for i = el_std(:)'; patch(mNdCrd(mLNodS(i,:),1),mNdCrd(mLNodS(i,:),2),'b','edgecolor','k','facecolor','y','facealpha',0.2); end
for i = el_hvi(:)'; patch(mNdCrd(mLNodS(i,:),1),mNdCrd(mLNodS(i,:),2),'b','edgecolor','k','facecolor','y','facealpha',1.0); end

nd_hvi = mLNodS(el_hvi,:); nd_hvi = unique(nd_hvi(mn_hvi));
nd_std = mLNodS(el_std,:); nd_std = unique(nd_std(mn_std));

scatter(mNdCrd(nd_std,1),mNdCrd(nd_std,2),50,'k')
scatter(mNdCrd(nd_hvi,1),mNdCrd(nd_hvi,2),40,'w')

%%