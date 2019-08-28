function SavePDF(name,h)

% Export a pdf figure without large margins

if nargin == 1 || isempty(h)
    h = gcf; % get current figure
end

% % get specific axes handle
% a = get(gcf,'Children');
% a = a(2); axes(a);

% get the current axes handle
a = get(h,'CurrentAxes');

if ~isempty(a)
    
    % get axes size
    set(a,'units','centimeters');
    pos = get(a,'OuterPosition');
    
    % resize paper so it fits axes
    set(h,'PaperUnits','centimeters');
    set(h,'PaperSize',[pos(3),pos(4)]);
    
    % change paper position
    set(h,'PaperPositionMode','manual');
    set(h,'PaperPosition',[0,0,pos(3),pos(4)]);
    
    % save
    saveas(h,[name,'.pdf']);
    
    % undo modifications to the figure
    set(a,'units','normalized');
    set(a,'OuterPosition',[0,0,1,1]);
    
    clear name_fig h a pos
    
else
    warning('Empty figure')
end