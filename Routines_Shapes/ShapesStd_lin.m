function [mGsShp,mGsDrv] = ShapesStd_lin(nGauss,mGsPnt)

nElNod = 2;

mGsShp(nGauss,nElNod) = 0;
mGsDrv(nGauss,nElNod) = 0;

for iGauss = 1:nGauss
    
    [mGsShp(iGauss,:),mGsDrv(iGauss,:)] = ...
        LgBasis_lin(mGsPnt(iGauss,:),nElNod);
    
end