function FigResize(szfig,szfnt,h)

% In case label moves outiside figure's outer position:
%   set(hx,'position',get(hx,'position')+[dx,dy,dz]);

if nargin == 2 || isempty(h)
    h = gcf; % get current figure
end

switch length(szfig)
    case 1
        szfig = [1,1,szfig,szfig];
    case 2
        szfig = [1,1,szfig];
    case 3
        szfig = [szfig,szfig(3)];
end

a = get(h,'CurrentAxes');
x = get(a,'XLabel');
y = get(a,'YLabel');
z = get(a,'ZLabel');
t = get(a,'Title');

set(a,'ActivePositionProperty','OuterPosition')

set(h,'Position',szfig);
set(a,'FontSize',szfnt);
set(x,'FontSize',szfnt);
set(y,'FontSize',szfnt);
set(z,'FontSize',szfnt);
set(t,'FontSize',szfnt);

end