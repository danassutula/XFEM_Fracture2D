function [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt)

nElNod = 4;
nDimes = 2;

mGsShp(nGauss,nElNod) = 0;
mGsDrv(nDimes*nGauss,nElNod) = 0;

jGsDrv = 1:nDimes;

for iGauss = 1:nGauss
    
    [mGsShp(iGauss,:),mGsDrv(jGsDrv,:)] = ...
        LgBasis_quad(mGsPnt(iGauss,:),nElNod);
    
    jGsDrv = jGsDrv + nDimes;
    
end