
%==========================================================================
% Plot Von Mises Stress (log-scale)
%==========================================================================


%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf;
set(h,'Color','w');

% dot size std.
dotsz_std = 6;
% dot size enr.
dotsz_enr = 1;

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

% vMises = [vMises_std;vMises_enr];
% mGsCrd = [mGsCrd_std;mGsCrd_enr];

%--------------------------------------------------------------------------
% Only plot what's within a bounding box (if at all specified)
%--------------------------------------------------------------------------
 
if exist('plotBox_stress','var') ...
        && ~isempty(plotBox_stress) ...
        && norm(diff(plotBox_stress,1))>0
    
    q_std = mGsCrd_std(:,1) > plotBox_stress(1,1) & mGsCrd_std(:,1) < plotBox_stress(2,1) & ...
            mGsCrd_std(:,2) > plotBox_stress(1,2) & mGsCrd_std(:,2) < plotBox_stress(2,2) ;
    
    q_enr = mGsCrd_enr(:,1) > plotBox_stress(1,1) & mGsCrd_enr(:,1) < plotBox_stress(2,1) & ...
            mGsCrd_enr(:,2) > plotBox_stress(1,2) & mGsCrd_enr(:,2) < plotBox_stress(2,2) ;
   
else
    
    q_std = 1:length(vMises_std);
    q_enr = 1:length(vMises_enr);
    
end

%--------------------------------------------------------------------------
% Plotting (with adjusted contrast)
%--------------------------------------------------------------------------

if exist('mGsDsp_std','var') && ...
        size(mGsCrd_std,1) == size(mGsDsp_std,1) && ...
        size(mGsCrd_enr,1) == size(mGsDsp_enr,1)
        
    mGsDfm_std = mGsCrd_std + mGsDsp_std*deformed_scale;
    mGsDfm_enr = mGsCrd_enr + mGsDsp_enr*deformed_scale;
    
    flag_plotCracks = 0;
    
else
    
    mGsDfm_std = mGsCrd_std;
    mGsDfm_enr = mGsCrd_enr;
    
    flag_plotCracks = 1;
    
end

clf(h); % better to clear the figure just before plotting
hold on;
axis equal;

scatter(mGsDfm_std(q_std,1),mGsDfm_std(q_std,2),...
    dotsz_std,log10(vMises_std(q_std)),'o','fill')

scatter(mGsDfm_enr(q_enr,1),mGsDfm_enr(q_enr,2),...
    dotsz_enr,log10(vMises_enr(q_enr)),'o','fill')

if flag_plotCracks && exist('cCkCrd','var')
    for ii = 1:length(cCkCrd)
        plot(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),'w','linewidth',1)
    end
end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

xmin = min([mGsDfm_std(q_std,1);mGsDfm_enr(q_enr,1)]);
xmax = max([mGsDfm_std(q_std,1);mGsDfm_enr(q_enr,1)]);
ymin = min([mGsDfm_std(q_std,2);mGsDfm_enr(q_enr,2)]);
ymax = max([mGsDfm_std(q_std,2);mGsDfm_enr(q_enr,2)]);

m = max(xmax-xmin,ymax-ymin)*0.05;
axis([xmin-m,xmax+m,ymin-m,ymax+m])
set(gca,'layer','top','box','on')

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
    title('von Mises stress')
end

% set color sinc.
colormap(jet) % bone, hot, cool
c = colorbar('location','EastOutside');
c.Label.String = ['log_{10} (\sigma_{vms})'];

% clear var.
clear vMises mGsCrd 
clear vMises_std vMises_enr
clear mGsDfm_std mGsDfm_enr 
clear xmin ymax ymin ymax m
clear flag_plotCracks 

figure(h);

%--------------------------------------------------------------------------
