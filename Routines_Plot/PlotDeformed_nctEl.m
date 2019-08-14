
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf; hold on;
set(h,'Color','w');

if ~exist('deformed_scale','var')
    deformed_scale = 10; end
if ~exist('deformed_dotsz','var')
    deformed_dotsz = 3; end

deformed_scale
deformed_dotsz

color_cut = [1.0,0.0,0.0];
color_nct = [0.9,0.9,0.9];
color_scl = sqrt(E(:)/max(E));

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

x_vrt = zeros(nLNodS,nElems); 
y_vrt = zeros(nLNodS,nElems);

for kk = 1:nElems
    x_vrt(:,kk) = mNdCrd(mLNodS(kk,:),1) + ...
        mNDspS(1,mLNodS(kk,:))' * deformed_scale;
    y_vrt(:,kk) = mNdCrd(mLNodS(kk,:),2) + ...
        mNDspS(2,mLNodS(kk,:))' * deformed_scale;
end

elems_cut = unique([cell2mat(cBrStd_eTp(:)); ...
    cell2mat(cHvBln_elm(:));cell2mat(cHvStd_elm(:))]);  % cut elements
elems_nct = setdiff([vElStd(:);vElEnr(:)],elems_cut);  % not cut elements

for ii = 1:nPhase
    
    q = intersect(elems_nct,find(vElPhz == ii))';
    patch(x_vrt(:,q),y_vrt(:,q),1+color_scl(ii,:)*(color_nct-1));
    
end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

axis equal

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
else
    xlabel('x')
    ylabel('y')
end

clear color_cut color_nct color_scl
clear elems_cut elems_nct
clear x_vrt y_vrt

%--------------------------------------------------------------------------
