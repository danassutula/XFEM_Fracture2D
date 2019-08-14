
%--------------------------------------------------------------------------
% Plot potential energy as a function of combined crack length
%--------------------------------------------------------------------------

% tip: do not clear the figure here in case you wish to superpose plots

h = gcf; hold on;
set(h,'Color','w');

% if iStep == nStep
%     % finished post-processing last solution increment
%     q = 1:iStep+1;
% else
%     % last material state was not post-processed; could be due:
%     % (1) a large energy change or (2) negative strain energy
%     q = 1:iStep;
% end

q = 1:min(length(vCkLen),length(Pi));

if ~isempty(q)
    h_plot = plot(vCkLen(q),Pi(q),'-sk', ...
        'linewidth',1,'markerfacecolor','r','markersize',5);
end

if length(q)>1
    set(gca,'XLim',[vCkLen(1),vCkLen(end)])
end

ylabel('Potential energy, \Pi')
xlabel('Combined crack surface, a')
title('Potential energy')

figure(h);

%--------------------------------------------------------------------------
