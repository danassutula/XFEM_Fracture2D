function [coord,lnods,nx] = StepDown_q4(xlim,ylim,nx,ny)

if mod(nx,2) == 0
    error('Number of nodes on the x-axis should be an odd number, i.e. mod(nx,2) == 1')
end

l = xlim(2)-xlim(1);
d = ylim(2)-ylim(1);

hx = l/(nx-1);         % initial size
hy = d/(2-1/2^(ny-2)); % initial size (by geometric progression)

lnods = [];
coord = [];

for j = 1:ny-1 % element row
    
    if j > 1
        ne = size(lnods,1); 
        nn = size(coord,1) - nx;
    else
        ne = 0;
        nn = 0;
    end
    
    lnods = [lnods;zeros(6*(nx-1)/2,4)];
    
for i = 1:(nx-1)/2 % element pairs (column wise)

    k = ne + (i-1)*6;
    
    % nd. numbering
    % 7--8--9-10-11
    % ---4--5--6---
    % 1-----2-----3
    
    n1 = nn + (i-1)*2+1;
    n2 = n1 + 1;
    n3 = n2 + 1;
    
    n4 = nn + nx + (i-1)*3+1;
    n5 = n4 + 1;
    n6 = n5 + 1;
    
    n7  = nn + nx + (nx-1)/2*3 + (i-1)*4+1;
    n8  = n7  + 1;
    n9  = n8  + 1;
    n10 = n9  + 1;
    n11 = n10 + 1;
    
    lnods(k+1,:) = [ n1, n2, n5, n4];
    lnods(k+2,:) = [ n4, n5, n9, n8];
    lnods(k+3,:) = [ n1, n4, n8, n7];
    lnods(k+4,:) = [ n2, n3, n6, n5];
    lnods(k+5,:) = [ n3,n11,n10, n6];
    lnods(k+6,:) = [ n5, n6,n10, n9];
    
end

xmid = [linspace(xlim(1)+hx/2  ,xlim(2)-hx*3/2,(nx-1)/2); ...
        linspace(xlim(1)+hx    ,xlim(2)-hx    ,(nx-1)/2); ...
        linspace(xlim(1)+hx*3/2,xlim(2)-hx/2  ,(nx-1)/2)];

if j > 1; ylim(1)=ylim(1)+2*hy; % advance the y-reference datum
    coord = [coord; [[xmid(:)',linspace(xlim(1),xlim(2),2*nx-1)]', ...
        [repmat(ylim(1)+hy/2,3*(nx-1)/2,1);repmat(ylim(1)+hy,2*nx-1,1)]]];
else
    coord = [[linspace(xlim(1),xlim(2),nx),xmid(:)',linspace(xlim(1),xlim(2),2*nx-1)]', ...
        [repmat(ylim(1),nx,1);repmat(ylim(1)+hy/2,3*(nx-1)/2,1);repmat(ylim(1)+hy,2*nx-1,1)]];
end

hx = hx/2;
hy = hy/2;

nx = 2*nx-1; % n. nodes on x-axis for the next layer of elements

end