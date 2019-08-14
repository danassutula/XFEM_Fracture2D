
%--------------------------------------------------------------------------
% PlotMesh
%--------------------------------------------------------------------------

h = gcf; clf(h);
set(h,'Color','w');
hold on; axis equal;

% color of stiffest mat.
c0 = [0.9,0.9,0.9];
% c0 = [1,1,1];

if exist('E','var')
    % sqrt-scaling color of all mat.
    c = 1+sqrt(E(:)/max(E))*(c0-1);
else
    c = repmat(c0,nElems,1);
end
 
if ~exist('nElems','var') || ~exist('nLNodS','var') 
    [nElems,nLNodS] = size(mLNodS);
end

x_vrt = zeros(nLNodS,nElems); 
y_vrt = zeros(nLNodS,nElems);

if exist('nPhase','var')
    for ii = 1:nPhase;
        q = find(vElPhz==ii)';
        for kk = q
            x_vrt(:,kk) = mNdCrd(mLNodS(kk,:),1);
            y_vrt(:,kk) = mNdCrd(mLNodS(kk,:),2);
        end
        patch(x_vrt(:,q),y_vrt(:,q),c(ii,:),'EdgeColor','k');
    end
else
    for kk = 1:nElems
        x_vrt(:,kk) = mNdCrd(mLNodS(kk,:),1);
        y_vrt(:,kk) = mNdCrd(mLNodS(kk,:),2);
    end
    patch(x_vrt,y_vrt,c(1,:),'EdgeColor','k');
end

% if nElems < 100
%     for ii = 1:nNdStd
%         text(mNdCrd(ii,1),mNdCrd(ii,2),num2str(ii),'FontSize',8,...
%             'color','r','BackgroundColor','w','edgecolor','r')
%     end
%     for ii = 1:nElems; q = mLNodS(ii,:);
%         text(sum(mNdCrd(q,1))/nLNodS,sum(mNdCrd(q,2))/nLNodS,num2str(ii),...
%             'FontSize',8,'color','b')
%     end
% end

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
