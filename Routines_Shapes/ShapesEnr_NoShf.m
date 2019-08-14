function [nLNodE,mGsShE,mGsDvE] = ShapesEnr_NoShf(...
    mElCrd,mGsShS,mGsDvS,nEnFun,mGFnVl,mGFnDv)
%==========================================================================
%
%
%==========================================================================

%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================
[nGauss,nLNodS] = size(mGsShS);
nLNodE = nLNodS*nEnFun;

nGsDrv = size(mGsDvS,1);
nShDrv = nGsDrv/nGauss;

mGsShE = zeros(nGauss,nLNodE);
mGsDvE = zeros(nGsDrv,nLNodE);

jGsDrv = 1:nShDrv;
jLNodS = 1:nLNodS;
jEnFun = 1:nEnFun;

for iGauss = 1:nGauss
    
    vElShp = mGsShS(iGauss,:);
    mElDrv = mGsDvS(jGsDrv,:);
    
    vFnVal = mGFnVl(iGauss,:);
    mFnDrv = mGFnDv(jGsDrv,:);
    
    mJackb = mElDrv*mElCrd;
    
    jLNodE = jLNodS;
    
    for iEnFun = jEnFun
        
        uFnVal = vFnVal(iEnFun);
        vFnDrv = mJackb*mFnDrv(:,iEnFun);
        
        mGsShE(iGauss,jLNodE) = vElShp.*uFnVal;
        mGsDvE(jGsDrv,jLNodE) = mElDrv.*uFnVal + vFnDrv*vElShp;
        
        jLNodE = jLNodE + nLNodS;
        
    end
    
    jGsDrv = jGsDrv + nShDrv;
    
end
%==========================================================================

end