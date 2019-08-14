function [vTpElm,vBlElm,mBlRmp] = Elems_tipBox(mNdCrd,mElNod,x_tip,r_tip)

%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================

mNdCrd(:,1) = mNdCrd(:,1) - x_tip(1);
mNdCrd(:,2) = mNdCrd(:,2) - x_tip(2);

mNdCrd = abs(mNdCrd*([1,1;-1,1]*0.5*sqrt(2)));
vRdNrm = mNdCrd(:,1)+mNdCrd(:,2);
vNdTip = find(vRdNrm < r_tip*(0.5*sqrt(2)));

mNdMbr = ismember(mElNod,vNdTip);
vNdSum = sum(mNdMbr,2);

vMbElm = find(vNdSum > 0);
vBlElm = vNdSum(vMbElm) < size(mElNod,2);

vTpElm = vMbElm(~vBlElm);
vBlElm = vMbElm( vBlElm);
mBlRmp = mNdMbr( vBlElm,:);

%==========================================================================

end