
%==========================================================================
% PlotDisplace
%==========================================================================


%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf;
set(h,'Color','w');

%--------------------------------------------------------------------------
% Get displacements
%--------------------------------------------------------------------------

mGsCrd = cell(nElems,1);
mGsDsp = cell(nElems,1);

for ii = vElStd(:)'
    
    vLNodS = mLNodS(ii,:);
    mLDspS = mNDspS(:,vLNodS);
    
    mGsCrd{ii} = mGsStd_omgShp*mNdCrd(vLNodS,:);
    mGsDsp{ii} = mGsStd_omgShp*mLDspS';
    
end

for ii = vElEnr(:)'
    
    vLNodS = mLNodS(ii,:);
    vLNodE = cLNodE{ii};
    
    mLDspS = mNDspS(:,vLNodS);
    mLDspE = mNDspE(:,vLNodE);
    
    mGsCrd{ii} = cGsEnr_omgShS{ii}*mNdCrd(vLNodS,:);
    mGsDsp{ii} = cGsEnr_omgShS{ii}*mLDspS'+cGsEnr_omgShE{ii}*mLDspE';
    
end

mGsCrd = cell2mat(mGsCrd);
mGsDsp = cell2mat(mGsDsp);

%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

clf(h); % better to clear the figure just before plotting
hold on;
axis equal;

% plot displacement vector field
% quiver(mGsCrd(:,1),mGsCrd(:,2),mGsDsp(:,1),mGsDsp(:,2),scale)

% scatter(mGsCrd(:,1),mGsCrd(:,1),3,mGsDsp(:,1),'fill')  
% scatter(mGsCrd(:,1),mGsCrd(:,2),3,mGsDsp(:,2),'fill')

scatter(mGsCrd(:,1),mGsCrd(:,2),50,...
    hypot(mGsDsp(:,1),mGsDsp(:,2)),'fill')

% plot cracks
for ii = 1:nCrack
    plot(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),'k','linewidth',1)
end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

xmin = min(mNdCrd(:,1));
xmax = max(mNdCrd(:,1));
ymin = min(mNdCrd(:,2));
ymax = max(mNdCrd(:,2));

m = max(xmax-xmin,ymax-ymin)*0.05;
axis([xmin-m,xmax+m,ymin-m,ymax+m])
set(gca,'layer','top','box','on')

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
end

colorbar

title('Displacement magnitude')

% clear var.
clear mGsCrd mGsDsp 
clear xmin ymax ymin ymax m

figure(h);

%--------------------------------------------------------------------------
