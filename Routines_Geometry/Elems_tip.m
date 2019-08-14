function [vTpElm,vBlElm,mBlRmp] = Elems_tip(mNdCrd,mElNod,x_tip,r_tip)

% Elems_tip:
%   returns elements, vTpElm and vBlElm, that have respectively all and
%   some of there nodes within a disk centered about x_tip of radius r_tip
%   mBlRmp is the blending ramp associated with elements vBlElm that is 1
%   in the interior of the disk and 0 on the otside; mNdCrd are the nodal 
%   coordinates; mElNod is the element topology.

mNdCrd(:,1) = mNdCrd(:,1) - x_tip(1);
mNdCrd(:,2) = mNdCrd(:,2) - x_tip(2);

vRdNrm = mNdCrd(:,1).^2+mNdCrd(:,2).^2;
vNdTip = find(vRdNrm < r_tip^2);

mNdMbr = ismember(mElNod,vNdTip);
vElMbr = find(any(mNdMbr,2));
vTpElm = all(mNdMbr(vElMbr,:),2);

vBlElm = vElMbr(~vTpElm);
mBlRmp = mNdMbr(vBlElm,:);
vTpElm = vElMbr( vTpElm);

end