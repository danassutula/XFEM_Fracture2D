
%--------------------------------------------------------------------------
% Plot energy dissipation rate: Gs = -dPi/da
%--------------------------------------------------------------------------

% tip: do not clear the figure here in case you wish to superpose plots

h = gcf; hold on;
set(h,'Color','w');

% if iStep == nStep
%     % finished post-processing last solution increment
%     q = 1:iStep;
% else
%     % last material state was not post-processed; could be due:
%     % (1) a large energy change or (2) negative strain energy
%     q = 1:iStep-1;
% end

q = 1:min(length(vCkLen)-1,length(Gs));

if ~isempty(q)
    h_plot = plot(0.5*(vCkLen(1:q(end))+vCkLen(2:q(end)+1)),Gs(q), ...
        '-^k','linewidth',1,'markerfacecolor','g','markersize',5);
end

if length(q)>1
    set(gca,'XLim',[vCkLen(1),vCkLen(q(end)+1)])
end

ylabel('Gs = -\Delta\Pi/\Deltaa')
xlabel('Combined crack surface, a')
title('Global energy dissipation')

figure(h);

%--------------------------------------------------------------------------
