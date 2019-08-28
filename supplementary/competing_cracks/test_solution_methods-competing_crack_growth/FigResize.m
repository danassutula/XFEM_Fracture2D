function FigResize(sz_fig,sz_fnt,h)

% In case label moves outiside figure's outer position:
%   set(hx,'position',get(hx,'position')+[dx,dy,dz]);

if nargin == 2 || isempty(h)
    h = gcf; % get current figure
end

switch length(sz_fig)
    case 1
        sz_fig = [1,1,sz_fig,sz_fig];
    case 2
        sz_fig = [1,1,sz_fig];
    case 3
        sz_fig = [sz_fig,sz_fig(3)];
end

a = get(h,'CurrentAxes');
x = get(a,'XLabel');
y = get(a,'YLabel');
z = get(a,'ZLabel');
t = get(a,'Title');

set(a,'ActivePositionProperty','OuterPosition')

set(h,'Position',sz_fig);
set(a,'FontSize',sz_fnt);
set(x,'FontSize',sz_fnt);
set(y,'FontSize',sz_fnt);
set(z,'FontSize',sz_fnt);
set(t,'FontSize',sz_fnt);

end