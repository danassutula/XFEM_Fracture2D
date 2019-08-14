function [LS_lin,LS_tp1,LS_tp2] = LevelSet2Line(x_pnt,x_lin)

%--------------------------------------------------------------------------
%
% OUTPUT:
%   LS_lin = signed distances normal to line segment ( + above, - bellow)
%   LS_tp1 = signed distances normal to line tip#1 ( + inside, - outside)
%
% INPUT:
%   x_pnt = matrix of point coordinates, e.g. rand(nPtCrd,2)
%   x_lin = line tip coordinates, e.g. rand(2,2)
%
%--------------------------------------------------------------------------

t1 = x_lin(2,1)-x_lin(1,1);
t2 = x_lin(2,2)-x_lin(1,2);

l = sqrt(t1^2+t2^2);

t1 = t1/l; n2 =  t1;
t2 = t2/l; n1 = -t2;

r1 = x_pnt(:,1)-x_lin(1,1);
r2 = x_pnt(:,2)-x_lin(1,2);

LS_lin = r1*n1+r2*n2;
LS_tp1 = r1*t1+r2*t2;
LS_tp2 = l-LS_tp1;

end