function [H,S,D,I] = LevelSet2Poly(x,p)

%--------------------------------------------------------------------------
% PART I (get distance)
%--------------------------------------------------------------------------

tol = 1e-13;

nx = size(x,1);
np = size(p,1);

jpA = 1:np-1;
jpB = jpA+1;

l1 = p(jpA,1)-p(jpB,1);
l2 = p(jpA,2)-p(jpB,2);

ll = l1.^2+l2.^2;

H(nx,1) = 0;
S(nx,1) = 0;
D(nx,2) = 0;
I(nx,1) = 0;

S = S + inf;

for ip = jpA
    
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
    jx = find(S-s>s*tol);
    
    D(jx,1) = d1(jx);
    D(jx,2) = d2(jx);
    
    S(jx) = s(jx);
    
    I(jx) = ip;
    
end

S = sqrt(S);

%--------------------------------------------------------------------------
% PART II (get sign)
%--------------------------------------------------------------------------

for ip = jpA

    jx = find(I==ip);
    
    H(jx) = D(jx,2)*l1(ip)-D(jx,1)*l2(ip);

end

jpA = 1:np-2;
jpB = jpA+1;

jpA = jpA(l1(jpA).*l1(jpB)+l2(jpA).*l2(jpB) < 0);     % get sharp kinks
jpB = jpA+1; crv = l2(jpA).*l1(jpB)-l1(jpA).*l2(jpB); % (-1)*curvature

for ic = 1:length(crv); ipA = jpA(ic); ipB = jpB(ic);
    
    jx = I==ipA & (p(ipB,1)-x(:,1))*l1(ipA)+(p(ipB,2)-x(:,2))*l2(ipA) > 0;
    
    H(jx) = crv(ic);
    
end

H(H>=0) = 1;
H(H< 0) =-1;

end
