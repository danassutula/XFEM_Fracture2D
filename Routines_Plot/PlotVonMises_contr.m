
%==========================================================================
% Plot Von Mises Stress (contours)
%==========================================================================


%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf; axis equal;
set(h,'Color','w');

% n. contours
ncont = 15;
% dot size
dotsz = 6;

if ~exist('deformed_scale','var')
    deformed_scale = 10;
end

%--------------------------------------------------------------------------
% Evaluate von Mises at Gauss points (if not done already)
%--------------------------------------------------------------------------

if ~exist('mStres_std','var') || isempty(mStres_std) ...
        || ~exist('mGsCrd_std','var') || isempty(mGsCrd_std)
    
    [mStres_std,mGsCrd_std] = Stress_std(mNDspS,mNdCrd,vElStd,mLNodS,...
        vElPhz,nPhase,cDMatx,mPrLod,mGsStd_omgShp,mGsStd_omgDrv);
    
    mStres_std = single(mStres_std);
    mGsCrd_std = single(mGsCrd_std);
    
    [mStres_enr,mGsCrd_enr] = Stress_enr(mNDspS,mNDspE,mNdCrd,vElEnr,mLNodS,cLNodE,...
        vElPhz,nPhase,cDMatx,mPrLod,cGsEnr_omgShS,cGsEnr_omgDvS,cGsEnr_omgDvE);

    mStres_enr = single(mStres_enr);
    mGsCrd_enr = single(mGsCrd_enr);
    
end

vMises_std = Stress_vms(mStres_std);
vMises_enr = Stress_vms(mStres_enr);

vMises = [vMises_std;vMises_enr];
mGsCrd = [mGsCrd_std;mGsCrd_enr];

%--------------------------------------------------------------------------
% Only plot what's within a bounding box (if at all specified)
%--------------------------------------------------------------------------

if exist('plotBox_stress','var') ...
        && ~isempty(plotBox_stress) ...
        && norm(diff(plotBox_stress,1))>0
    
    q = mGsCrd(:,1) > plotBox_stress(1,1) & mGsCrd(:,1) < plotBox_stress(2,1) & ...
        mGsCrd(:,2) > plotBox_stress(1,2) & mGsCrd(:,2) < plotBox_stress(2,2) ;
    
    vMises = vMises(q);
    mGsCrd = mGsCrd(q,:);
    
else
    
    q = 1:length(vMises);
    
end

%--------------------------------------------------------------------------
% Adjust contrast
%--------------------------------------------------------------------------

% % set sig. bands [min, max]
% sig_bnd = linspace(min(vMises),max(vMises),ncont+1);
% sig_bnd = 0.5*(sig_bnd(1:end-1)+sig_bnd(2:end));

% set sig. band: [0, 1]
sig_bnd = linspace(0,1,ncont+1);
sig_bnd = 0.5*(sig_bnd(1:end-1)+sig_bnd(2:end));

[~,ii] = sort(vMises);

n = length(ii);

if ncont > n
    ncont = n;
end
    
dn = round(n/ncont-0.5);
kk  = 1:(n-dn*(ncont-1));
vMises(ii(kk)) = sig_bnd(1);

kk = length(kk)+(1:dn);
for jj = 2:ncont
    vMises(ii(kk)) = sig_bnd(jj); kk = kk + dn;
end

%--------------------------------------------------------------------------
% Plotting (with adjusted contrast)
%--------------------------------------------------------------------------

if exist('mGsDsp_std','var') && ...
        size(mGsCrd_std,1) == size(mGsDsp_std,1) && ...
        size(mGsCrd_enr,1) == size(mGsDsp_enr,1)
    
    mGsDsp = [mGsDsp_std;mGsDsp_enr]*deformed_scale;
    mGsDfm = mGsCrd(q,:) + mGsDsp(q,:);
    
    flag_plotCracks = 0;
    
else
    
    mGsDfm = mGsCrd(q,:);
    
    flag_plotCracks = 1;
    
end

clf(h); % better to clear the figure just before plotting
hold on;
axis equal;

scatter(mGsDfm(:,1),mGsDfm(:,2),dotsz,vMises,'o','fill')

if flag_plotCracks && exist('cCkCrd','var')
    for ii = 1:length(cCkCrd)
        plot(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),'w','linewidth',1)
    end
 end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

xmin = min(mGsDfm(:,1));
xmax = max(mGsDfm(:,1));
ymin = min(mGsDfm(:,2));
ymax = max(mGsDfm(:,2));

m = max(xmax-xmin,ymax-ymin)*0.05;
axis([xmin-m,xmax+m,ymin-m,ymax+m])
set(gca,'layer','top','box','on')

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
    title('von Mises stress contours')
end

% set color sinc.
colormap(jet) % bone, hot, cool
colorbar('location','EastOutside')

% clear var.
clear vMises mGsCrd mGsDfm
clear vMises_std vMises_enr
clear xmin ymax ymin ymax m
clear flag_plotCracks 

figure(h);

%--------------------------------------------------------------------------
