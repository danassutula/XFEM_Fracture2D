function [mGsPnt,vGsWgt] = Gauss_DomQuad(nGauss)

%==========================================================================
% Linear domain is transformation to Quadratic domain
%==========================================================================

nGauss = sqrt(nGauss);

[vGsPtL,vGsWgL] = Gauss_DomLin(nGauss);

mGsPtX = vGsPtL(:,ones(1,nGauss));
mGsPtY = mGsPtX';

mGsPnt = [mGsPtX(:),mGsPtY(:)];

mGsWgt = vGsWgL*vGsWgL';
vGsWgt = mGsWgt(:);

end

%==========================================================================
