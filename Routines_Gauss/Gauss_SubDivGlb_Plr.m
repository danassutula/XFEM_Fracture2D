function [nGauss_mod,mGsPnt_mod,vGsWgt_mod] = ...
    Gauss_SubDivGlb_Plr(mElCrd,mNdCrd,mSbNod,nGsPlr)
%==========================================================================

%--------------------------------------------------------------------------
% OUTPUT:
%
% mGsPnt_mod - modified Gauss points
% vGsWgt_mod - modified weights
%
% INPUT:
%
% mElCrd - element coordinates (geometric typically)
% mNdCrd - global coordinates of sub-elements
% mSbNod - topology of sub-elements
% nGsPlr - [nGsTht,nGsRdi] number of Gauss points
%--------------------------------------------------------------------------

%==========================================================================
%                           BULLETPROOF
%==========================================================================

%--------------------------------------------------------------------------
% PART I
%--------------------------------------------------------------------------

[nGauss_mod,mGsPnt_mod,vGsWgt_mod] = ...
    Gauss_SubDivLcl_Plr(mNdCrd,mSbNod,nGsPlr);

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