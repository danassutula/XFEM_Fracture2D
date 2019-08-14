
%--------------------------------------------------------------------------
% PlotDomain
%--------------------------------------------------------------------------

h = gcf; clf(h);
set(h,'Color','w');
hold on; axis equal;

% color of stiffest material
c0 = [0.8,0.8,0.8];

% sqrt-scaling color of all mat.
c = 1+sqrt(E(:)/max(E))*(c0-1);

x_vrt = zeros(nLNodS,nElems); 
y_vrt = zeros(nLNodS,nElems);

for ii = 1:nPhase; q = find(vElPhz==ii)';
    
    for kk = q
        x_vrt(:,kk) = mNdCrd(mLNodS(kk,:),1);
        y_vrt(:,kk) = mNdCrd(mLNodS(kk,:),2);
    end
    
    patch(x_vrt(:,q),y_vrt(:,q),c(ii,:),'EdgeColor',c(ii,:));
    
end

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

% clear var.'s
clear c0 c x_vrt y_vrt 
clear xmin ymax ymin ymax m

figure(h);

%--------------------------------------------------------------------------
