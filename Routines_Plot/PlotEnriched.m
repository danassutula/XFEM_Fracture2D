
%==========================================================================
% Plot Enriched Elements
%==========================================================================

h = gcf; % clf(h);
set(h,'Color','w');
hold on; axis equal;

facealpha = 0.3;

x_vrt = zeros(nLNodS,nElEnr);
y_vrt = zeros(nLNodS,nElEnr);

%--------------------------------------------------------------------------
% Branch (std.) - 1st layer
%--------------------------------------------------------------------------

n = 0;

for ii = 1:nCrack
    for jj = 1:2
        for kk = [cBrStd_eTp{ii,jj};cBrStd_eSp{ii,jj};cBrStd_eFl{ii,jj}]'
            
            n = n + 1;
            
            x_vrt(:,n) = mNdCrd(mLNodS(kk,:),1);
            y_vrt(:,n) = mNdCrd(mLNodS(kk,:),2);
            
        end
    end
end

patch(x_vrt(:,1:n),y_vrt(:,1:n),[0,0,1],'facealpha',facealpha);

%--------------------------------------------------------------------------
% Branch (bln.) - 2nd layer
%--------------------------------------------------------------------------

n = 0;

for ii = 1:nCrack
    for jj = 1:2
        for kk = [cBrBln_eSp{ii,jj};cBrBln_eFl{ii,jj}]'
            
            n = n + 1;
            
            x_vrt(:,n) = mNdCrd(mLNodS(kk,:),1);
            y_vrt(:,n) = mNdCrd(mLNodS(kk,:),2);
            
        end
    end
end

patch(x_vrt(:,1:n),y_vrt(:,1:n),[0,1,1],'facealpha',facealpha);

%--------------------------------------------------------------------------
% Heaviside (all) - 3st layer
%--------------------------------------------------------------------------

n = 0;

for ii = 1:nCrack
    for jj = cHvStd_elm{ii}'
        
        n = n + 1;
        
        x_vrt(:,n) = mNdCrd(mLNodS(jj,:),1);
        y_vrt(:,n) = mNdCrd(mLNodS(jj,:),2);
        
    end
end

patch(x_vrt(:,1:n),y_vrt(:,1:n),[1,0,0],'facealpha',facealpha); % red

n = 0;

for ii = 1:nCrack
    for jj = [cHvBln_elm{ii,1};cHvBln_elm{ii,2}]'
        
        n = n + 1;
        
        x_vrt(:,n) = mNdCrd(mLNodS(jj,:),1);
        y_vrt(:,n) = mNdCrd(mLNodS(jj,:),2);
        
    end
end

patch(x_vrt(:,1:n),y_vrt(:,1:n),[1,1,0],'facealpha',facealpha); % yellow

%--------------------------------------------------------------------------
% Plot SIF evaluation radius
%--------------------------------------------------------------------------

for ii = 1:nCrack
    PlotCircle(mTpRdi(ii,1)*f_sif,cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'--y');
    PlotCircle(mTpRdi(ii,2)*f_sif,cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'--y');
end

%--------------------------------------------------------------------------
% Plot crack tip intersection criterion radius
%--------------------------------------------------------------------------

for ii = 1:nCrack
    PlotCircle(mTpRdi(ii,1)*f_xrs,cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'-.r');
    PlotCircle(mTpRdi(ii,2)*f_xrs,cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'-.r');
end

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

% clear var.'s
clear facealpha
clear x_vrt y_vrt
clear xmin ymax ymin ymax m

figure(h);
