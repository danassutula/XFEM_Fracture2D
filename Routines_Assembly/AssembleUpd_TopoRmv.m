
%==========================================================================
% Remove former enrichment data
%==========================================================================

% real to logical
if ~islogical(vElEnr); p = vElEnr;
    vElEnr = false(nElems,1);
    vElEnr(p) = true;
end

% % n. nd. enr. (upper bnd.)
% nNdEn8 = size(mNDofE,2);

% init. el. & nd. to rmv.
vElRmv = false(nElems,1);
vNdRmv = false(nNdEn8,1);

%--------------------------------------------------------------------------
% Branch enrichment
%--------------------------------------------------------------------------

for i = 1:2 % tip #i
for j = find(mCk2Up(:,i))' % crack #j
    
    % Brn. el. to rmv.
    elrmv = [ ...
    cBrStd_eTp{j,i};
    cBrStd_eSp{j,i};
    cBrStd_eFl{j,i};
    cBrBln_eFl{j,i};
    cBrBln_eSp{j,i} ];
    
    % rmv. ck. intersection (EXCEPTION TO DO IT HERE!)
    for k = cBrStd_eTp{j,i}'
        cXsElm{k} = [];
    end
    
    % clr. Brn. el.
    cBrStd_eTp{j,i} = [];
    cBrStd_eSp{j,i} = [];
    cBrStd_eFl{j,i} = [];
    cBrBln_eSp{j,i} = [];
    cBrBln_rSp{j,i} = [];
    cBrBln_eFl{j,i} = [];
    cBrBln_rFl{j,i} = [];
    
    % el. to rmv.
    vElRmv(elrmv) = true;
    % nd. to rmv.
    vNdRmv(cNdEnr_brn{j,i}) = true;
    % clr. Brn. nd.
    cNdEnr_brn{j,i} = [];
    
    for k = elrmv(:)' % rmv. shapes and nodes only for this enr.
        
        % el. enr. data
        mLEnDt = cLEnDt{k};
        vLNodE = cLNodE{k};
        
        % get enrichment position in mLEnDt
        p0 = find(mLEnDt(1,:)==j & mLEnDt(2,:)==i);
        
        % remove corresponding enrichment data
        cLEnDt{k}(:,p0) = [];
        
        % get position outside Brn. enr.
        p = sum(mLEnDt(3,1:p0-1));
        p = [1:p,p+mLEnDt(3,p0)+1:length(vLNodE)];
        
        % enr. nodes to keep
        cLNodE{k} = vLNodE(p);
        
        % enr. shapes to keep
        cGsEnr_omgDvE{k} = cGsEnr_omgDvE{k}(:,p);
        cGsEnr_omgShE{k} = cGsEnr_omgShE{k}(:,p);
        
    end
    
end
end

%--------------------------------------------------------------------------
% Some post-processing
%--------------------------------------------------------------------------

% mod. enr. el.
vElRmv = find(vElRmv);
nElRmv = length(vElRmv);

% org. el. to upd.
vElUp0 = false(nElRmv,1);

for i = 1:nElRmv % cond.: any remaining enr.
    vElUp0(i) = size(cLEnDt{vElRmv(i)},2) > 0;
end

% current enr. el.
vElEnr(vElRmv(~vElUp0)) = false;

% el. to upd.
vElUp0 = vElRmv(vElUp0);
nElUp0 = length(vElUp0);

%--------------------------------------------------------------------------
% Remove former DOF and resize Kg, Fg
%--------------------------------------------------------------------------

% enr. nd. to rmv.
vNdRmv = find(vNdRmv);

% get dofs to remove
p = mNDofE(:,vNdRmv);
p = p(:);

% rmv. dofs from Kg
mGlStf(:,p) = [];
mGlStf(p,:) = [];
% rmv. dofs from Fg
vGlFrc(p,:) = [];

% set unused dofs to zero
mNDofE(:,vNdRmv) = 0;

% get currently enr. nodes
vNdEnr = find(mNDofE(1,:));
nNdEnr = length(vNdEnr);

% get current n. of dof.
nGDofE = nDimes*nNdEnr;
nGlDof = nGDofS+nGDofE;

for i = 1:nDimes % shift dofs since some were removed
    mNDofE(i,vNdEnr) = nGDofS+i:nDimes:nGlDof-nDimes+i;
end

%--------------------------------------------------------------------------
% Updated Kg, Fg
%--------------------------------------------------------------------------

if nElUp0 > 0; p = vElUp0;
    
    % subtract Ke (se,es,ee)
    mGlStf = mGlStf - MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),...
    mNDofE,vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));

    % % subtract Fe (ck. face)
    % if wCkLod
    %   'Fg is to be recomputed from scratch'
    % end
    
    % subtract Fe (pre-load)
    if wPrLod
        vGlFrc = vGlFrc - ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),...
        cLNodE(p),mNDofE,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));
    end
    
    % subtract Fe (body)
    if wBdLod
        vGlFrc = vGlFrc - ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),...
        cLNodE(p),mNDofE,cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);
    end
    
end

%--------------------------------------------------------------------------
% Remove shapes and weights
%--------------------------------------------------------------------------

% does not work this way:
% cGsEnr_omgShS(vElRmv) = [];

for i = vElRmv(:)' % contains vElUp0
    
    % remove shapes (std.)
    cGsEnr_omgShS{i} = [];
    cGsEnr_omgDvS{i} = [];
    % remove shapes (enr.)
    cGsEnr_omgDvE{i} = [];
    cGsEnr_omgShE{i} = [];
    % remove weights
    cGsEnr_omgWgt{i} = [];
    
end

%--------------------------------------------------------------------------
% Delete temporary variables
%--------------------------------------------------------------------------

clear elrmv

%--------------------------------------------------------------------------
