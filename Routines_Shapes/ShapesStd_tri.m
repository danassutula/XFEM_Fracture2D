function [mGsShp,mGsDrv] = ShapesStd_tri(nGauss,mGsPnt)

nElNod = 3;
nDimes = 2;

mGsShp(nGauss,nElNod) = 0;
mGsDrv(nDimes*nGauss,nElNod) = 0;

jGsDrv = 1:nDimes;

for iGauss = 1:nGauss
    
    [mGsShp(iGauss,:),mGsDrv(jGsDrv,:)] = ...
        LgBasis_tri(mGsPnt(iGauss,:),nElNod);
    
    jGsDrv = jGsDrv + nDimes;
    
end