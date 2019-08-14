function [coord,lnods,nx] = StepUp_t3(xlim,ylim,nx,ny)

l = xlim(2)-xlim(1);
d = ylim(2)-ylim(1);

hx = l/(nx-1);         % initial size
hy = d/(2^(ny-1)-1)/2; % initial size (by geometric progression)

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
    
    lnods = [lnods;zeros(3*(nx-1)/2,3)];
    
for i = 1:(nx-1)/2 % element pairs (column wise)
    
    k = ne + 3*(i-1);
    
    % nd. numbering
    % 4_____5
    % |\   /|
    % 1_\2/_3
    
    n1 = nn + (i-1)*2+1;
    n2 = n1 + 1;
    n3 = n2 + 1;
    
    n4 = nn + nx + i;
    n5 = n4 + 1;
    
    lnods(k+1,:) = [ n1, n2, n4];
    lnods(k+2,:) = [ n2, n5, n4];
    lnods(k+3,:) = [ n2, n3, n5];
    
end

if j > 1; ylim(1)=ylim(1)+hy; % advance the y-reference datum
    coord = [coord; [[linspace(xlim(1),xlim(2),(nx-1)/2+1)]', ...
        [repmat(ylim(1)+2*hy,(nx-1)/2+1,1)]]];
else
    coord = [[linspace(xlim(1),xlim(2),nx),linspace(xlim(1),xlim(2),(nx-1)/2+1)]', ...
        [repmat(ylim(1),nx,1);repmat(ylim(1)+2*hy,(nx-1)/2+1,1)]];
end

hx = hx*2;
hy = hy*2;

nx = (nx-1)/2+1; % n. nodes on x-axis for the next layer of elements

end