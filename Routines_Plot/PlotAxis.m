function axbox = PlotAxis(margin,x)
%==========================================================================
%
% Get position of axsis
%
%==========================================================================

xmrg = (max(max(x)) - min(min(x)))*margin; % 10% margin
xbar = mean(x);

x(:,1) = x(:,1) - xbar(1);
x(:,2) = x(:,2) - xbar(2);

xmin = xbar(1) + min(x(:,1)) - xmrg;
xmax = xbar(1) + max(x(:,1)) + xmrg;
ymin = xbar(2) + min(x(:,2)) - xmrg;
ymax = xbar(2) + max(x(:,2)) + xmrg;

axbox = [xmin,xmax,ymin,ymax];

%==========================================================================

end