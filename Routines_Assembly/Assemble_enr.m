
%==========================================================================
% Assemble global (enr.)
%==========================================================================


%--------------------------------------------------------------------------
% Enrichment control
%--------------------------------------------------------------------------

nEnFun_hvi = 1;
nEnFun_brn = 3;

nLNodE_brn = nLNodS*nEnFun_brn;
nLNodE_hvi = nLNodS*nEnFun_hvi;

nEnJmp_brn = 1;
jEnJmp_brn = 1; % 1st enr. fun has jump (but, in general, can be a vector)

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Determine enriched element topology
%--------------------------------------------------------------------------
AssembleEnr_Topology
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Determine standard shape functions for enriched elements
%--------------------------------------------------------------------------
AssembleEnr_ShapesStd
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Determine enriched shape functions for enriched elements
%--------------------------------------------------------------------------
AssembleEnr_ShapesEnr
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Ensure size compatibility of enr. shape functions and enr. nodes
%--------------------------------------------------------------------------

for i = vElEnr(:)'
    if size(cLNodE{i},2) ~= size(cGsEnr_omgShE{i},2)
        error(['Missmatch of size compatibility of ',...
            'enr. shape functions and enr. nodes. Element #',num2str(i)])
    end
end

%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Get enr. dofs
%--------------------------------------------------------------------------

vNdEnr = 1:nNdEnr;

nGDofE = nNdEnr*nDimes;
nGlDof = nGDofS+nGDofE;

mNDofE = zeros(nDimes,nNdEnr);
mNDofE(:) = nGDofS+1:nGlDof;

if with_Update % saves mem.
    mNDofE = sparse(mNDofE);
end

%--------------------------------------------------------------------------
% Assemble Kg (enr.)
%--------------------------------------------------------------------------

mGlStf = MtxStiffEnr(nGlDof,mNdCrd,mLNodS,vElEnr,cLNodE,mNDofE,...
    vElPhz,nPhase,cDMatx,cGsEnr_omgDvS,cGsEnr_omgDvE,cGsEnr_omgWgt);

%--------------------------------------------------------------------------
% Assemble Fg (enr.)
%--------------------------------------------------------------------------

vGlFrc = sparse([],[],[],nGlDof,1);

if wCkLod
    
    vGlFrc_crk = ...
        ... % Heaviside Elements (optimized)
        ForceEnr_BndCkHvi(nGlDof,nCrack,cCkCrd, ...
        mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod) +   ...
        ... % Branch Elements (not optimized)
        ForceEnr_BndCkBrn(nGlDof,nElEnr,jEnJmp_brn,nLNodS,nCrack,cCkCrd,...
        nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,mCkLod);
        
    % for now keep crack face force separate
    vGlFrc = vGlFrc + vGlFrc_crk;
    
end

if wPrLod
    vGlFrc = vGlFrc + ForceEnr_Resid(nGlDof,mNdCrd,mLNodS,vElEnr,...
        cLNodE,mNDofE,cGsEnr_omgDvS,cGsEnr_omgDvE,cGsEnr_omgWgt,mPrLod);
end

if wBdLod
    vGlFrc = vGlFrc + ForceEnr_Body(nGlDof,mNdCrd,mLNodS,vElEnr,...
        cLNodE,mNDofE,cGsEnr_omgDvS,cGsEnr_omgShE,cGsEnr_omgWgt,vBdLod);
end
