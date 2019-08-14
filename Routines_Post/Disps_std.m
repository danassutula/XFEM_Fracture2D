function [mGsDsp_std,mGsCrd_std] = Disps_std(vElStd)

%--------------------------------------------------------------------------
% Declare Global Variables
%--------------------------------------------------------------------------

global mNDspS
global mNdCrd
global mLNodS
global mGsStd_omgShp

%--------------------------------------------------------------------------

nElStd = length(vElStd);
mGsCrd_std = cell(nElStd,1);
mGsDsp_std = cell(nElStd,1);

for i = 1:nElStd
    
    vLNodS = mLNodS(vElStd(i),:);
    
    mGsCrd_std{i} = mGsStd_omgShp*mNdCrd(vLNodS,:);
    mGsDsp_std{i} = mGsStd_omgShp*mNDspS(:,vLNodS)';
    
end

mGsCrd_std = cell2mat(mGsCrd_std);
mGsDsp_std = cell2mat(mGsDsp_std);