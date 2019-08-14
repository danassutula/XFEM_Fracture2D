function [S,D,I] = LevelSet2Poly_abs(x,p)

tol = 1e-14;
tol = tol^2;

nx = size(x,1);
np = size(p,1);

jp = 1:np-1;

l1 = p(jp,1)-p(jp+1,1);
l2 = p(jp,2)-p(jp+1,2);

ll = l1.^2+l2.^2;

S(nx,1) = 0;
D(nx,2) = 0;
I(nx,1) = 0;

S = S + inf;

for ip = jp
    
    d1 = p(ip,1)-x(:,1);
    d2 = p(ip,2)-x(:,2);
    
    ld = d1*l1(ip)+d2*l2(ip);
    
    jx = find(ld < ll(ip) & ld > 0);
    
    f = ld(jx)./ll(ip);
    
    d1(jx) = d1(jx)-l1(ip)*f;
    d2(jx) = d2(jx)-l2(ip)*f;
    
    jx = find(ld >= ll(ip));
    
    d1(jx) = d1(jx)-l1(ip);
    d2(jx) = d2(jx)-l2(ip);
    
    s = d1.^2+d2.^2;
    
    % comparing squares
    jx = find(S-s > tol);
    
    D(jx,1) = d1(jx);
    D(jx,2) = d2(jx);
    
    S(jx) = s(jx);
    
    I(jx) = ip;
    
end

S = sqrt(S);