
%--------------------------------------------------------------------------
% Plot Cracks (iterations)
%--------------------------------------------------------------------------

h = gcf; hold on;
set(h,'Color','w');
% axis equal;

for ii = find(mCk2Up(:,1))'
    
    plot(cCkCrd{ii}([1,2],1),cCkCrd{ii}([1,2],2),':or', ...
        'linewidth',0.8,'markersize',3,'markerfacecolor','r')
    
end

for ii = find(mCk2Up(:,2))'
    
    plot(cCkCrd{ii}([end,end-1],1),cCkCrd{ii}([end,end-1],2),':or', ...
        'linewidth',0.8,'markersize',3,'markerfacecolor','r')
    
end

figure(h);

%--------------------------------------------------------------------------
