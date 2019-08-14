
%==========================================================================
% Shapes for Branch enriched elements (bln.) 
%==========================================================================

for iElBrn = 1:length(vElBrn)
    
    uElBrn = vElBrn(iElBrn);
    vBlRmp = mBlRmp(iElBrn,:);
    
    mElCrd = mNdCrd(mLNodS(uElBrn,:),:);
    
    %----------------------------------------------------------------------
    % Gauss shapes (std.)
    %----------------------------------------------------------------------
    
    mGsShS = cGsEnr_omgShS{uElBrn};
    mGsDvS = cGsEnr_omgDvS{uElBrn};
    
    %----------------------------------------------------------------------
    % Gauss shapes (enr.)
    %----------------------------------------------------------------------
    
    mGsCrd = mGsShS*mElCrd;
    
    [mGFnVl,mGFnDv,mNdShf] = EnrFun_Brn(...
        mElCrd,mGsCrd,nEnFun_brn,vTpCrd,uTpAlf);
    
    vGBlVl = mGsShS*vBlRmp(:);
    vGBlDv = mGsDvS*vBlRmp(:);
    
    [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtBln(mElCrd,mGsShS,...
        mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf,vGBlVl,vGBlDv);
    
    cGsEnr_omgShE{uElBrn}(:,end+1:end+nLNodE) = mGsShE;
    cGsEnr_omgDvE{uElBrn}(:,end+1:end+nLNodE) = mGsDvE;
    
    %----------------------------------------------------------------------
    
end