function [mGsCrd,vGsWgt] = Gauss_DomTri_Plr(mVxCrd,nGsPlr)
%==========================================================================
%
%
%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================

nGsTht = nGsPlr(1);
nGsRdi = nGsPlr(2);

%--------------------------------------------------------------------------
[vGsPtT,vGsWtT] = Gauss_DomLin(nGsTht); % vGsPtT = (-1..1);
[vGsPtR,vGsWtR] = Gauss_DomLin(nGsRdi); % vGsPtR = (-1..1);

vThShp = 0.5*(1+vGsPtT(:));
vRdShp = 0.5*(1+vGsPtR(:));
%--------------------------------------------------------------------------

nGauss = nGsTht*nGsRdi;

mGsCrd(nGauss,2) = 0;
vGsRdi(nGauss,1) = 0;
mRdDir(nGsTht,2) = 0;

v1 = mVxCrd(2,:)-mVxCrd(1,:); u1 = sqrt(v1(1)^2+v1(2)^2)\v1;
v2 = mVxCrd(3,:)-mVxCrd(1,:); u2 = sqrt(v2(1)^2+v2(2)^2)\v2;

dv = v2 - v1;

uThLim = acos(u1(1)*u2(1)+u1(2)*u2(2)); % angle between the two vectors
vGsTht = vThShp*uThLim+atan2(v1(2),v1(1)); % OR GetPol(v1(1),v1(2));

mRdDir(:,1) = cos(vGsTht);
mRdDir(:,2) = sin(vGsTht);

iGauss = 1:nGsRdi;
for iGsTht = 1:nGsTht
    
    r_dir = mRdDir(iGsTht,:);
    r_lim = (dv(1)*v1(2)-dv(2)*v1(1))/(dv(1)*r_dir(2)-dv(2)*r_dir(1));
    
    vGsJac = vRdShp*r_lim; % cartesian/polar jackobian
    
    vGsRdi(iGauss)   = vGsJac*r_lim; % times by r_lim so to avoid doing later 
    mGsCrd(iGauss,:) = vGsJac*r_dir;
    
    iGauss = iGauss + nGsRdi;
    
end

mGsCrd(:,1) = mGsCrd(:,1) + mVxCrd(1,1);
mGsCrd(:,2) = mGsCrd(:,2) + mVxCrd(1,2);

dJ = 0.25*uThLim; % 0.5*0.5 -> come from |J_tht| and |J_rdi| 
vGsWgt = vGsRdi.*kron(vGsWtT,vGsWtR)*dJ; % Kronecker tensor product

%==========================================================================

end