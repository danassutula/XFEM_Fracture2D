function [nCrack,cCkCrd] = GenCrack_layer(nCrack,a,xstd,xlim,ylim,tlim)
%
% nCrack = number of cracks
% cCkCrd = cell of crack coordinates
%
% a     = length of every crack
% xstd  = std. deviation in crack spacing
% xlim  = distribute cracks on 'x' interval
% ylim  = distribute cracks within 'y' interval
% tlim  = interval for crack rotation angle

lt = xlim(2)-xlim(1);
lc = a*nCrack;
lg = lt-lc;

if lg < 0
    warning('Negative crack spacing; the maximum number of cracks, n_max = %d',nCrack-ceil(abs(lg)/a));
end

d = lg/nCrack;

x = linspace(xlim(1)+0.5*(d+a),xlim(2)-0.5*(d+a),nCrack) + xstd*randn(1,nCrack);
y = ylim(1) + (ylim(2)-ylim(1))*rand(1,nCrack);
t = tlim(1) + (tlim(2)-tlim(1))*rand(1,nCrack);

dx = 0.5*cos(t)*a;
dy = 0.5*sin(t)*a;

cCkCrd = cell(nCrack,1);

for i = 1:nCrack
    
    x2 = x(i) + dx(i); 
    y2 = y(i) + dy(i);
    
    x1 = x(i) - dx(i); 
    y1 = y(i) - dy(i);
    
    cCkCrd{i} = [x1,y1;x2,y2];

end
