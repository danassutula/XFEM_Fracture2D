
%--------------------------------------------------------------------------
% Plot Cracks
%--------------------------------------------------------------------------

h = gcf; hold on;
set(h,'Color','w');
axis equal;

param_PlotCracks_markerType  = 'o';
param_PlotCracks_lineColor   = 'k';
param_PlotCracks_markerColor = 'w';

% param_PlotCracks_markerType  = '+';
% param_PlotCracks_lineColor   = 'r';
% param_PlotCracks_markerColor = 'r';

%--------------------------------------------------------------------------
% Plot Cracks
%--------------------------------------------------------------------------

% handle for each plot
h_plot = zeros(nCrack,1);

for ii = 1:nCrack
    
    h_plot(ii) = plot(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),...
        ['-',param_PlotCracks_markerType,param_PlotCracks_lineColor], ...
        'markerfacecolor',param_PlotCracks_markerColor,...
        'markeredgecolor',param_PlotCracks_lineColor,...
        'linewidth',1,'markersize',5);
    
end

%--------------------------------------------------------------------------
% Mark crack tip activity
%--------------------------------------------------------------------------

if exist('mTpAct','var') && size(mTpAct,1) == nCrack
    
    % active crack tips (green)
    for ii = find(mTpAct(:,1))'
        plot(cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'o','markersize',5,...
            'markerfacecolor','g','markeredgecolor','k')
    end
    for ii = find(mTpAct(:,2))'
        plot(cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'o','markersize',5,...
            'markerfacecolor','g','markeredgecolor','k')
    end
    
    % inactive crack tips (red)
    if exist('mCkJun','var') && size(mCkJun,1) == nCrack
        for ii = find(~mTpAct(:,1) & ~mCkJun(:,1))'
            plot(cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'o','markersize',5,...
                'markerfacecolor','r','markeredgecolor','k')
        end
        for ii = find(~mTpAct(:,2) & ~mCkJun(:,2))'
            plot(cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'or','markersize',5,...
                'markerfacecolor','r','markeredgecolor','k')
        end
    end
    
end

%--------------------------------------------------------------------------

if exist('mNdCrd','var')
    
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
    clear x y xmin ymax ymin ymax m
    
end

figure(h);

%--------------------------------------------------------------------------
