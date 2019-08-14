function l = ElemAnticlock(l,x)
%--------------------------------------------------------------------------
%
% ElemAnticlock: 
%   Make element node order anti-clockwise (positive Jackobian)
%
%--------------------------------------------------------------------------

[ne,nn] = size(l);

qL = 2:nn;
qR = nn:-1:2;

for i = 1:ne
    if (x(l(i,2),1)-x(l(i,1),1))*(x(l(i,3),2)-x(l(i,1),2)) - ...
            (x(l(i,3),1)-x(l(i,1),1))*(x(l(i,2),2)-x(l(i,1),2)) < 0
       
       l(i,qL) = l(i,qR);
       
    end
end