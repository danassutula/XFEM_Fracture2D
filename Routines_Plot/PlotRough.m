
%--------------------------------------------------------------------------
% PlotRough
%--------------------------------------------------------------------------

h = gcf; clf(h);
set(h,'Color','w');

plot(xi,yi,'k',[xi(1),xi(end)],[y0,y0],':k')

axis normal
set(gca,'box','on');

xlabel('Position, x')
ylabel('Profile, y')

legend(num2str(Rq,'y, (R_q = %0.3e)'),'y_{mean}','location','SouthEast');

%--------------------------------------------------------------------------
