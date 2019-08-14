function [p,t,b] = MeshRect_t3(xlim,ylim,h0)
    
L = xlim(2)-xlim(1);
H = ylim(2)-ylim(1);

nx = round(L/h0);
ny = round(H/(h0*sqrt(3)/2));

if ny-2*floor(ny/2) == 1
   ny = ny+1; % make even 
end

np = (nx+1)*(ny+1)+ny/2;
ne = (2*nx+1)*ny/2;

p = zeros(np,2);
t = zeros(ne,3);
b = cell(4,1);

hx = L/nx;
hy = H/ny;

x1 = xlim(1):  hx:xlim(2);
y1 = ylim(1):2*hy:ylim(2);

x2 =[xlim(1),x1(1:end-1)+hx/2,xlim(2)];
y2 = y1(1:end-1)+hy;

Tp = 2*(nx+1)+1; % N points per layer
Te = 2*(2*nx+1); % N elements per layer

for i = 1:ny/2 % go by period
    
    % node coordinates for this layer
    
    p((i-1)*Tp+1:(i-1)*Tp+nx+1,1) = x1;
    p((i-1)*Tp+1:(i-1)*Tp+nx+1,2) = y1(i);
    
    p(i*Tp-nx-1:i*Tp,1) = x2;
    p(i*Tp-nx-1:i*Tp,2) = y2(i);
    
    % element topology for this layer (right-angle triangles)
    
    t((i-1)*Te+1,:)     = [(i-1)*Tp+1    ,i*Tp-nx ,i*Tp-nx-1]; % SW
    t((i-0.5)*Te,:)     = [    i*Tp-nx-2 ,i*Tp    ,i*Tp-1   ]; % SE
    
    t((i-0.5)*Te+1,:)   = [    i*Tp-nx-1 ,i*Tp-nx ,i*Tp+1   ]; % NW
    t(i*Te,:)           = [    i*Tp-1    ,i*Tp    ,i*Tp+nx+1]; % NE
    
    % element topology for this layer (equilateral triangles)
    
    q = (i-1)*Te+2:2:(i-0.5)*Te-1;
    
    t(q,1) = t((i-1)*Te+1,1)    : t((i-0.5)*Te,1) -1; % SS
    t(q,2) = t((i-1)*Te+1,1) +1 : t((i-0.5)*Te,1)   ; % SS
    t(q,3) = t((i-1)*Te+1,2)    : t((i-0.5)*Te,3)   ; % SS
    
    q = (i-1)*Te+3:2:(i-0.5)*Te-2;
    
    t(q,1) = t((i-1)*Te+1,1) +1 : t((i-0.5)*Te,1) -1; % S
    t(q,2) = t((i-1)*Te+1,2) +1 : t((i-0.5)*Te,3)   ; % S
    t(q,3) = t((i-1)*Te+1,2)    : t((i-0.5)*Te,3) -1; % S
    
    q = (i-0.5)*Te+3:2:i*Te-2;
    
    t(q,1) = t((i-0.5)*Te+1,2)    : t(i*Te,1) -1; % N
    t(q,2) = t((i-0.5)*Te+1,2) +1 : t(i*Te,1)   ; % N
    t(q,3) = t((i-0.5)*Te+1,3) +1 : t(i*Te,3) -1; % N
    
    q = (i-0.5)*Te+2:2:i*Te-1;
    
    t(q,1) = t((i-0.5)*Te+1,2)    : t(i*Te,1)   ; % NN
    t(q,2) = t((i-0.5)*Te+1,3) +1 : t(i*Te,3)   ; % NN
    t(q,3) = t((i-0.5)*Te+1,3)    : t(i*Te,3) -1; % NN
    
    
end

% node coordinates for top boundary

p(np-nx:np,1) = x1;
p(np-nx:np,2) = y1(end);

% boundary nodes

b{1} = (1:nx+1)';
b{3} = (np-nx:np)';

b{2} = find(p(:,1) > xlim(2)-h0/4);
b{4} = find(p(:,1) < xlim(1)+h0/4);
