
%==========================================================================
% Shapes for Branch enriched elements 
%==========================================================================

for uElBrn = vElBrn(:)'
    
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
    
    [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtShf(mElCrd,...
        mGsShS,mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf);
    
    cGsEnr_omgShE{uElBrn}(:,end+1:end+nLNodE) = mGsShE;
    cGsEnr_omgDvE{uElBrn}(:,end+1:end+nLNodE) = mGsDvE;
    
    %----------------------------------------------------------------------
    
end