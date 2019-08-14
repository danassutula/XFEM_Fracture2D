function [nGauss_mod,mGsPnt_mod,vGsWgt_mod] = ...
    Gauss_SubDivLcl_Plr(mNdCrd,mElNod,nGsPlr)
%==========================================================================

%--------------------------------------------------------------------------
% OUTPUT:
%
% mGsPnt_mod - modified Gauss points
% vGsWgt_mod - modified weights
%
% INPUT:
%
% mNdCrd - local coordinates of sub-elements
% mElNod - topology of sub-elements
% nGsPlr - [nGsTht,nGsRdi] number of Gauss points
%--------------------------------------------------------------------------

%==========================================================================
%                           BULLETPROOF
%==========================================================================

nElems_sub = size(mElNod,1);
nGauss_sub = nGsPlr(1)*nGsPlr(2);
nGauss_mod = nGauss_sub*nElems_sub;

mGsPnt_mod(nGauss_mod,2) = 0;
vGsWgt_mod(nGauss_mod,1) = 0;

iGauss_mod = 1:nGauss_sub;
for iElemn = 1:nElems_sub
    
    [mGsPnt_mod(iGauss_mod,:),vGsWgt_mod(iGauss_mod)] = ...
        Gauss_DomTri_Plr(mNdCrd(mElNod(iElemn,:),:),nGsPlr);
    
    iGauss_mod = iGauss_mod + nGauss_sub;
    
end

%==========================================================================
end