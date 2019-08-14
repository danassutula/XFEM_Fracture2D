function [vElMbr,mElWgt] = Elems_sif(mNdCrd,mElNod,x_tip,R)

mNdCrd(:,1) = mNdCrd(:,1) - x_tip(1);
mNdCrd(:,2) = mNdCrd(:,2) - x_tip(2);

vRadii = mNdCrd(:,1).^2+mNdCrd(:,2).^2;
mElWgt = ismember(mElNod,find(vRadii < R^2));

vElMbr = find(any(mElWgt,2) & ~all(mElWgt,2));
mElWgt = mElWgt(vElMbr,:)'; % (need transpose)

end