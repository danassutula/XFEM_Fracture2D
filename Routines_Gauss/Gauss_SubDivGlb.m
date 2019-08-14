function [nGauss_mod,mGsPnt_mod,vGsWgt_mod] = Gauss_SubDivGlb(...
    mElCrd,mNdCrd,mSbNod,nGauss_sub,mGsShp_sub,mGsDrv_sub,vGsWgt_sub)
%==========================================================================

%--------------------------------------------------------------------------
% PART I
%--------------------------------------------------------------------------

[nGauss_mod,mGsPnt_mod,vGsWgt_mod] = Gauss_SubDivLcl(...
    mNdCrd,mSbNod,nGauss_sub,mGsShp_sub,mGsDrv_sub,vGsWgt_sub);

%--------------------------------------------------------------------------
% PART II: map global Gauss points to local
%--------------------------------------------------------------------------

mGsPnt_mod = Gauss_Glb2Lcl(mGsPnt_mod,mElCrd);

%--------------------------------------------------------------------------
% PART III: modify weights
%--------------------------------------------------------------------------

nElNod = size(mElCrd,1);

for i = 1:nGauss_mod
    
    [~,dNdX] = LgBasis_omg(mGsPnt_mod(i,:),nElNod); dxdX = dNdX*mElCrd;
    vGsWgt_mod(i) = vGsWgt_mod(i)/(dxdX(1)*dxdX(4)-dxdX(2)*dxdX(3));
    
end

%==========================================================================
end