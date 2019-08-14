function [coord,lnods,bound] = MeshRect_q4(xlim,ylim,nx,ny)

nnode = nx*ny;
nelem = (nx-1)*(ny-1);

coord = zeros(nnode,2);
lnods = zeros(nelem,4);
bound = cell(4,1);

x = linspace(xlim(1),xlim(2),nx);
x = x(ones(ny,1),:)'; x = x(:);

y = linspace(ylim(1),ylim(2),ny);
y = y(ones(nx,1),:); y = y(:);

coord(:,1) = x;
coord(:,2) = y;

ie = 1;
for iy = 1:(ny-1)
    for ix = 1:(nx-1);
       
        lnods(ie,1) = (iy-1)*nx+ix;
        lnods(ie,2) = lnods(ie,1)+1;
        lnods(ie,3) = lnods(ie,2)+nx;
        lnods(ie,4) = lnods(ie,3)-1;
        
        ie=ie+1;
        
    end
end

bound{1} = (1:nx)';
bound{2} = ((1:ny)*nx)';
bound{3} = (ny-1)*nx+bound{1};
bound{4} = bound{2}-(nx-1);
