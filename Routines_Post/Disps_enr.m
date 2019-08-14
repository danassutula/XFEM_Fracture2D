function [mGsDsp_enr,mGsCrd_enr] = Disps_enr(vElEnr)

%--------------------------------------------------------------------------
% Declare Global Variables
%--------------------------------------------------------------------------

global mNDspS
global mNDspE
global mNdCrd
global mLNodS
global cLNodE
global cGsEnr_omgShS
global cGsEnr_omgShE

%--------------------------------------------------------------------------

nElEnr = length(vElEnr);
mGsCrd_enr = cell(nElEnr,1);
mGsDsp_enr = cell(nElEnr,1);

for i = 1:nElEnr; ii=vElEnr(i);
    
    if ~(isempty(cGsEnr_omgShS{ii}) || isempty(cGsEnr_omgShE{ii}))
        mGsCrd_enr{i} = cGsEnr_omgShS{ii}*mNdCrd(mLNodS(ii,:),:);
        mGsDsp_enr{i} = cGsEnr_omgShS{ii}*mNDspS(:,mLNodS(ii,:))' ...
            + cGsEnr_omgShE{ii}*mNDspE(:,cLNodE{ii})';
    else
        warning('enriched element''s shape functions are empty.')
    end 
    
end

mGsCrd_enr = cell2mat(mGsCrd_enr);
mGsDsp_enr = cell2mat(mGsDsp_enr);

