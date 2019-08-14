
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
vNdRmv = false(nNdEn8,1); % (to remove nodal dofs from eq.)

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
    
    vElRmv(elrmv) = true; % el. to rmv.
    vNdRmv(cNdEnr_brn{j,i}) = true; % nd. to rmv.
    cNdEnr_brn{j,i} = []; % clr. Brn. nd.
    
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
% Heaviside enrichment
%--------------------------------------------------------------------------

% ck. opposite end
aux = [2,1];

% el. eq. already rmv.
elup0 = false(nElems,1);

for i = 1:2
for j = find(mCk2Up(:,i))'
    
    % end-segment elements and end-kink elements will be removed
    
    if i == 1
        
        % tip & prev. sgm.
        x_inc = cCkCrd{j}([3,2,1],:);
        p_inc = cHvStd_sgm{j}(:,1) == 1;
        
    else
        
        % tip & prev. sgm.
        x_inc = cCkCrd{j}(end-2:end,:);
        p_inc = cHvStd_sgm{j}(:,2) == max(cHvStd_sgm{j}(:,2));
        
    end
    
    elhld = [cHvStd_elm{j}(~p_inc); cHvBln_elm{j,aux(i)}];
    elrmv = [cHvStd_elm{j}( p_inc); cHvBln_elm{j,i}     ];
    
    % rmv. el. from ck. inc.
    cHvStd_elm{j}(p_inc)   = [];
    cHvStd_sgm{j}(p_inc,:) = [];
    
    % rmv. bln. el.
    cHvBln_elm{j,i} = [];
    cHvBln_sgm{j,i} = [];
    cHvBln_rmp{j,i} = [];
    
    % el. to rmv./recomp.
    vElRmv(elrmv) = true;
    
    % init. Hvi. nd. to rmv.
    ndrmv = false(nNdEn8,1);
    
    % ck. inc. el. node members for finding nd. to rmv.
    ndmbr = ismember(mLNodS(elrmv,:),mLNodS(elhld,:));
    
    % DO NOT DO elvtx NOW (done subsequentlly)
    for kk = 1:length(elrmv); k = elrmv(kk);
        
        % el. enr. data
        mLEnDt = cLEnDt{k};
        vLNodE = cLNodE{k};
        
        % get enrichment position in mLEnDt
        p0 = find(mLEnDt(1,:)==j & mLEnDt(2,:)==0);
        
        % remove corresponding enrichment data
        cLEnDt{k}(:,p0) = [];
        
        % get position inside Hvi. enr.
        p = sum(mLEnDt(3,1:p0-1));
        p = p+1:p+mLEnDt(3,p0);
        
        % decouple elrmv from elhld
        if any(ndmbr(kk,:)) && ~elup0(k);
            
            % el. contrib. rmv.
            elup0(k) = true;
            
            % Hvi. nd. to remove: a subset
            ndrmv(vLNodE(p(~ndmbr(kk,:)))) = true;
            
            % subtract Ke (se,es,ee)
            mGlStf = mGlStf - MtxStiffEnr(nGlDof,mNdCrd,mLNodS(k,:),1,cLNodE(k),mNDofE,...
            vElPhz(k),nPhase,cDMatx,cGsEnr_omgDvS(k),cGsEnr_omgDvE(k),cGsEnr_omgWgt(k));
        
            % % subtract Fe (ck. face)
            % if wCkLod
            %   'Fg is to be recomputed from scratch'
            % end
            
            % subtract Fe (pre-load)
            if wPrLod
                vGlFrc = vGlFrc - ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(k,:),1,cLNodE(k),...
                mNDofE,cGsEnr_omgDvS(k),cGsEnr_omgDvE(k),cGsEnr_omgWgt(k),mPrLod(:,k));
            end
            
            % subtract Fe (body)
            if wBdLod
                vGlFrc = vGlFrc - ForceEnr_Body(nGlDof,mNdCrd,mLNodS(k,:),1,cLNodE(k),...
                mNDofE,cGsEnr_omgDvS(k),cGsEnr_omgShE(k),cGsEnr_omgWgt(k),vBdLod);
            end
            
        else
            
            % Hvi. nd. to remove
            ndrmv(vLNodE(p)) = true;
            
        end
        
        % get position outside Hvi. enr.
        p = [1:p(1)-1,p(end)+1:length(vLNodE)];
        
        % enr. nodes to keep
        cLNodE{k} = vLNodE(p);
        
        % enr. shapes to keep
        cGsEnr_omgDvE{k} = cGsEnr_omgDvE{k}(:,p);
        cGsEnr_omgShE{k} = cGsEnr_omgShE{k}(:,p);
        
        % Hvi. nd. to remove (glb.)
        vNdRmv(ndrmv) = true;
        
        % rmv. Hvi. nd. for this ck.
        p =~ismember(cNdEnr_hvi{j},find(ndrmv));
        cNdEnr_hvi{j} = cNdEnr_hvi{j}(p);
        
        % rmv. ck. xrs.
        cXsElm{k} = [];
        
    end
    
    % typially overwritten later, but try not to carry garbage
    cNdBln_hvi{j,i} = []; % (yet, will exist in cNdEnr_hvi{j,i})
    
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

% el. whose contirb. still need to rmv.
p = vElUp0(~ismember(vElUp0,find(elup0)));

if ~isempty(p)
    
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

clear aux elup0 elhld elrmv ndmbr ndrmv

%--------------------------------------------------------------------------
