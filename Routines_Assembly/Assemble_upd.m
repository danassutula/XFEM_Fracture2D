
% key variables used for updating: 
%   mCk2Up, mNoInc_grw, mNoInc_upd

% variables not used:
%   vTp2Up

% The following variables will not be changed:
%   mCk2Up nTp2Up vTp2Up mNoInc_grw

%--------------------------------------------------------------------------
% Update Fg before anything else (quick & dirty)
%--------------------------------------------------------------------------

if wCkLod % subtract ck. face lod.
    vGlFrc = vGlFrc - vGlFrc_crk;
end

%--------------------------------------------------------------------------


%==========================================================================
% Assemble global (upd.)
%==========================================================================

% n. nd. enr. (upper bnd.)
nNdEn8 = size(mNDofE,2);

% for tracking
nNdEn0 = nNdEn8;
nNdE00 = nNdEn8;

%--------------------------------------------------------------------------
% Backup/restore time-state
%--------------------------------------------------------------------------

if fIter_inc && ~iIter_dir
   if iIter_inc % ~= 0
       AssembleUpd_Restore
   else
       AssembleUpd_Backup
   end
end

%--------------------------------------------------------------------------
% Remove former enrichment and determine elements to update
%--------------------------------------------------------------------------

if iIter_dir % > 0 || == -1
    AssembleUpd_TopoRmv_itr
else % if iIter_dir == 0
    AssembleUpd_TopoRmv
end

%--------------------------------------------------------------------------
% Update crack segment numbering
%--------------------------------------------------------------------------

for i = find(mCk2Up(:,1) & mNoInc_grw(:,1))'
    
    cHvBln_sgm{i,1} = cHvBln_sgm{i,1} + mNoInc_grw(i,1); % push forward
    cHvStd_sgm{i}   = cHvStd_sgm{i}   + mNoInc_grw(i,1);
    cHvBln_sgm{i,2} = cHvBln_sgm{i,2} + mNoInc_grw(i,1);
    
end

%--------------------------------------------------------------------------
% Refine mesh if necessary
%--------------------------------------------------------------------------

% set crack tips to refine
mCk2Rf = mCk2Up & mCkRfn;

if any(mCk2Rf(:))
    
    AssembleUpd_TopoRfn
    
    % (Hvi. enrichment is done)
    mNoInc_upd(mCk2Rf) = 0;
    
    % no further refinements if fully coarsened
    mCkRfn(mCkRfN==0) = false;
    
end

%--------------------------------------------------------------------------
% Add new enrichment and determine elements to update
%--------------------------------------------------------------------------

AssembleUpd_TopoAdd

%--------------------------------------------------------------------------
% Update shapes (std.)
%--------------------------------------------------------------------------

AssembleUpd_ShapesStd

%--------------------------------------------------------------------------
% Update shapes (enr.)
%--------------------------------------------------------------------------

AssembleUpd_ShapesEnr

%--------------------------------------------------------------------------
% Update enr. dofs
%--------------------------------------------------------------------------

% get upd. enr. nd.
vNdUpd = nNdE00+1:nNdEn8;
nNdUpd = length(vNdUpd);

% upd. enr. nd.
nNdEnr = nNdEnr+nNdUpd;
vNdEnr = [vNdEnr,vNdUpd];

% upd. n. of enr. dofs
nGDofE = nNdEnr*nDimes;

% upd. tot. n. of dofs
nGlDf0 = nGlDof;
nGlDof = nGDofS+nGDofE;

% upd. dofs mtx.
for i = 1:nDimes
    mNDofE(i,vNdUpd) = nGlDf0+i:nDimes:nGlDof-nDimes+i;
end

% resize Kg & Fg
if nGlDof > nGlDf0
    mGlStf(nGlDof,nGlDof) = 0;
    vGlFrc(nGlDof) = 0;
end

%--------------------------------------------------------------------------


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Safety checks
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if with_Debug
    
    elbad_ch1 = [];
    elbad_ch2 = [];
    elbad_ch3 = [];
    elbad_ch4 = [];
    elbad_ch5 = [];
    elbad_ch6 = [];
    
    for i = vElStd(:)'
        if ~(isempty(cLNodE{i}) && isempty(cLEnDt{i}) && isempty(cXsElm{i}))
            elbad_ch1 = [elbad_ch1;i];
        end
    end
    
    for i = vElEnr(:)'
        
        if isempty(cLNodE{i}) || isempty(cLEnDt{i})
            elbad_ch2 = [elbad_ch2;i];
        end
        
        if length(cLNodE{i}) ~= sum(cLEnDt{i}(3,:))
            elbad_ch3 = [elbad_ch3;i];
        end
        
        if size(unique(cLEnDt{i}(1:2,:)','rows'),1) ~= size(cLEnDt{i},2)
            elbad_ch4 = [elbad_ch4;i];
        end
        
        if any(reshape(~mNDofE(:,cLNodE{i}),[],1))
            elbad_ch5 = [elbad_ch5;i];
        end
        
        if size(cLNodE{i},2) ~= size(cGsEnr_omgShE{i},2)
            elbad_ch6 = [elbad_ch6;i];
        end
        
    end
    
    ckbad_ch1 = [];
    
    for i = find(mTpAct(:,1))'
       if isempty(cHvBln_elm{i,1}) 
           ckbad_ch1 = [ckbad_ch1;i,1];
       end
    end
    for i = find(mTpAct(:,2))'
       if isempty(cHvBln_elm{i,2}) 
           ckbad_ch1 = [ckbad_ch1;i,2];
       end
    end
    
    fail = 0;
    
    if ~isempty(elbad_ch1); fail = 1;
        warning('(ch1)')
    end
    if ~isempty(elbad_ch2); fail = 1;
        warning('(ch2)')
    end
    if ~isempty(elbad_ch3); fail = 1;
        warning('(ch3)')
    end
    if ~isempty(elbad_ch4); fail = 1;
        warning('(ch4)')
    end
    if ~isempty(elbad_ch5); fail = 1;
        warning('(ch5)')
    end
    if ~isempty(elbad_ch6); fail = 1;
        warning('missmatch of size compatibility of enr. shape functions and enr. nodes (ch6)')
    end
    
    if ~isempty(ckbad_ch1); fail = 1;
        warning('no Heaviside blending elements were detected for the crack tips that are branch enriched')
        % go to script AssebleUpd_TopoRfn and uncomment lines 1081-1090 to
        % enable the correction (re-comment lines 1093-1098)
    end
    
    if fail
        error('safety check(s) failed')
    end

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if nElUpd > 0; p = vElUpd;
%--------------------------------------------------------------------------
% Update Kg (enr.)
%--------------------------------------------------------------------------
    
mGlStf = mGlStf + ...
    MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
    vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));

%--------------------------------------------------------------------------
% Update Fg (enr.)
%--------------------------------------------------------------------------

% if wCkLod
%     'force due to crack face tractions is recomputed'
% end

if wPrLod

    vGlFrc = vGlFrc + ...
        ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
        cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));

end

if wBdLod

    vGlFrc = vGlFrc + ...
        ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
        cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);

end
end

%--------------------------------------------------------------------------
% Update Fg (add crack surface tractions)
%--------------------------------------------------------------------------

if wCkLod
    
    vGlFrc_crk = ...
        ... % Heaviside Elements (optimized)
        ForceEnr_BndCkHvi(nGlDof,nCrack,cCkCrd, ...
        mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod) +   ...
        ... % Branch Elements (not optimized)
        ForceEnr_BndCkBrn(nGlDof,nElEnr,jEnJmp_brn,nLNodS,nCrack,cCkCrd,...
        nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,mCkLod);
    
    vGlFrc = vGlFrc + vGlFrc_crk;
    
end

%--------------------------------------------------------------------------
