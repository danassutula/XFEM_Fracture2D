function [coord,lnods,nx] = StepDown_t3(xlim,ylim,nx,ny)

% step down element size

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
    
    lnods = [lnods;zeros(3*(nx-1),3)];
    
for i = 1:nx-1 % element pairs (column wise)

    k = ne + 3*(i-1); % element caunter
    
    % nd. numbering
    
    % 3__4__5
    % | / \ |
    % 1/___\2
    
    n1 = nn + i;
    n2 = n1 + 1;
    
    n3 = nn + nx + (i-1)*2+1;
    n4 = n3 + 1;
    n5 = n4 + 1;
    
    lnods(k+1,:) = [ n1, n4, n3];
    lnods(k+2,:) = [ n1, n2, n4];
    lnods(k+3,:) = [ n2, n5, n4];
    
end

if j > 1; ylim(1)=ylim(1)+2*hy; % advance the y-reference datum
    coord = [coord; [linspace(xlim(1),xlim(2),2*nx-1)', ...
        [repmat(ylim(1)+hy,2*nx-1,1)]]];
else
    coord = [[linspace(xlim(1),xlim(2),nx),linspace(xlim(1),xlim(2),2*nx-1)]', ...
        [repmat(ylim(1),nx,1);repmat(ylim(1)+hy,2*nx-1,1)]];
end

hx = hx/2;
hy = hy/2;

nx = 2*nx-1; % n. nodes on x-axis for the next layer of elements

end