
%==========================================================================
% Plot Enriched Elements
%==========================================================================

h = gcf; hold on
set(gcf,'color','w')

%--------------------------------------------------------------------------
% Plot SIF evaluation radius
%--------------------------------------------------------------------------

for ii = 1:nCrack
    PlotCircle(mTpRdi(ii,1)*f_sif,cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'-k');
    PlotCircle(mTpRdi(ii,2)*f_sif,cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'-k');
end

% for ii = 1:nCrack
%     PlotCircle(mTpRdi(ii,1)*f_rmv,cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),'-.r');
%     PlotCircle(mTpRdi(ii,2)*f_rmv,cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),'-.r');
% end

%--------------------------------------------------------------------------
% Plot tip activity
%--------------------------------------------------------------------------

for ii = find(mTpAct(:,1))'
    scatter(cCkCrd{ii}(1,1),cCkCrd{ii}(1,2),20,'r','fill','markeredgecolor','k')
end

for ii = find(mTpAct(:,2))'
    scatter(cCkCrd{ii}(end,1),cCkCrd{ii}(end,2),20,'b','fill','markeredgecolor','k')
end

%--------------------------------------------------------------------------
% Plot tip inactivity
%--------------------------------------------------------------------------
% 
% for ii = find(~mTpAct(:,1))'
%     scatter(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),200,'r','fill')
% end
% 
% for ii = find(~mTpAct(:,1))'
%     scatter(cCkCrd{ii}(:,1),cCkCrd{ii}(:,2),200,'r','fill')
% end
%
%--------------------------------------------------------------------------