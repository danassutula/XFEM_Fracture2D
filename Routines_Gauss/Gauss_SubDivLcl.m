function [nGauss_mod,mGsPnt_mod,vGsWgt_mod] = Gauss_SubDivLcl(...
    mNdCrd,mElNod,nGauss_sub,mGsShp_sub,mGsDrv_sub,vGsWgt_sub)
%==========================================================================

nElems_sub = size(mElNod,1);
nGauss_mod = nGauss_sub*nElems_sub;

mGsPnt_mod(nGauss_mod,2) = 0;
vGsWgt_mod(nGauss_mod,1) = 0;

j_sub = 1:nGauss_sub;
j_sd0 = 1:2;
j_mod = 1;

for i = 1:nElems_sub
    
    x_elm = mNdCrd(mElNod(i,:),:);
    dxdX  = mGsDrv_sub*x_elm;
    
    j_sbd = j_sd0;
    for j = j_sub
        
        detJ = dxdX(j_sbd(1),1)*dxdX(j_sbd(2),2) - ...
               dxdX(j_sbd(2),1)*dxdX(j_sbd(1),2);
        
        mGsPnt_mod(j_mod,:) = mGsShp_sub(j,:)*x_elm;
        vGsWgt_mod(j_mod)   = vGsWgt_sub(j)*detJ;
        
        j_sbd = j_sbd + 2;
        j_mod = j_mod + 1;
        
    end
    
end

%==========================================================================
end