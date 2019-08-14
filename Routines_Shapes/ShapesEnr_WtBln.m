function [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtBln(mElCrd,...
    mGsShS,mGsDvS,nEnFun,mGFnVl,mGFnDv,mFnShf,vGBlVl,vGBlDv)
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

jFnRpt = ones(1,nShDrv);

if isempty(mFnShf)
    mFnShf = zeros(nEnFun,nLNodS);
end

switch isempty(mGFnDv)
    case 0
        
        mJackb = mGsDvS*mElCrd;
        
        for iGauss = 1:nGauss
            
            vElShp = mGsShS(iGauss,:);
            mElDrv = mGsDvS(jGsDrv,:);
            
            vFnVal = mGFnVl(iGauss,:);
            mFnDrv = mGFnDv(jGsDrv,:);
            
            uBlVal = vGBlVl(iGauss);
            vBlDrv = vGBlDv(jGsDrv);
            
            jLNodE = jLNodS;
            
            for iEnFun = jEnFun
                
                vFnShf = vFnVal(iEnFun) - mFnShf(iEnFun,:);
                vFnDrv = mJackb(jGsDrv,:)*mFnDrv(:,iEnFun);
                
                mGsShE(iGauss,jLNodE) = uBlVal*vElShp.*vFnShf;
                mGsDvE(jGsDrv,jLNodE) = mElDrv.*vFnShf(jFnRpt,:)*uBlVal + vFnDrv*uBlVal*vElShp + vBlDrv*(vElShp.*vFnShf);
                
                jLNodE = jLNodE + nLNodS;
                
            end
            
            jGsDrv = jGsDrv + nShDrv;
            
        end
        
    otherwise % case true
        
        for iGauss = 1:nGauss
            
            vElShp = mGsShS(iGauss,:);
            mElDrv = mGsDvS(jGsDrv,:);
            
            vFnVal = mGFnVl(iGauss,:);
            
            uBlVal = vGBlVl(iGauss);
            vBlDrv = vGBlDv(jGsDrv);
            
            jLNodE = jLNodS;
            
            for iEnFun = jEnFun
                
                vFnShf = vFnVal(iEnFun) - mFnShf(iEnFun,:);
                
                mGsShE(iGauss,jLNodE) = uBlVal*vElShp.*vFnShf;
                mGsDvE(jGsDrv,jLNodE) = mElDrv.*vFnShf(jFnRpt,:)*uBlVal + vBlDrv*(vElShp.*vFnShf);
                
                jLNodE = jLNodE + nLNodS;
                
            end
            
            jGsDrv = jGsDrv + nShDrv;
            
        end
        
end
%==========================================================================

end