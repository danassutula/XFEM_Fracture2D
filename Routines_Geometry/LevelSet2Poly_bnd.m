function [H,S,D,I] = LevelSet2Poly_bnd(x,p)

% get level set to a boundary that is defined by a closed polyline (first 
% and last nodes correspond); otherwise get level set to a polyline that is
% not closed

%--------------------------------------------------------------------------
% PART I: get distance
%--------------------------------------------------------------------------

tol = 1e-13;
fac = 1+tol;

nx = size(x,1);
np = size(p,1);

qA = 1:np-1;
qB = qA+1;

l1 = p(qA,1)-p(qB,1);
l2 = p(qA,2)-p(qB,2);

ll = l1.^2+l2.^2;

H(nx,1) = 0;
S(nx,1) = 0;
D(nx,2) = 0;
I(nx,1) = 0;

S = S + inf;

for i = qA
    
    d1 = p(i,1)-x(:,1);
    d2 = p(i,2)-x(:,2);
    
    ld = d1*l1(i)+d2*l2(i);
    
    qx = ld < ll(i) & ld > 0;
    
    f = ld(qx)./ll(i);
    
    d1(qx) = d1(qx)-l1(i)*f;
    d2(qx) = d2(qx)-l2(i)*f;
    
    qx = ld >= ll(i);
    
    d1(qx) = d1(qx)-l1(i);
    d2(qx) = d2(qx)-l2(i);
    
    s = d1.^2+d2.^2;
    
    % comp. squares
    qx = S > s*fac;
    
    D(qx,1) = d1(qx);
    D(qx,2) = d2(qx);
    
    S(qx) = s(qx);
    
    I(qx) = i;
    
end

S = sqrt(S);

%--------------------------------------------------------------------------
% PART II: get sign
%--------------------------------------------------------------------------

for i = qA; qx = I == i;    
    H(qx) = D(qx,2)*l1(i)-D(qx,1)*l2(i);
end

qA = 1:np-2;
qB = qA+1;

qA = qA(l1(qA).*l1(qB)+l2(qA).*l2(qB) < 0); % get sharp kinks
qB = qA+1; K=l2(qA).*l1(qB)-l1(qA).*l2(qB); % (-1)*curvature

for i = 1:length(K)
    
    jA = qA(i);
    jB = qB(i);
    
    qx = I==jA & (p(jB,1)-x(:,1))*l1(jA)+(p(jB,2)-x(:,2))*l2(jA)>0;
    
    H(qx) = K(i);
    
end

if p(1,:) == p(np,:) % (closed boundary)
    
    jA = np-1;
    jB = 1;
    
    if l1(jA).*l1(jB)+l2(jA).*l2(jB) < 0 % if sharp kink
            
        qx = I==jB; % 2nd 'qx' indicates association with p(1,:)
        qx = qx((p(jB,1)-x(qx,1))*l1(jB)+(p(jB,2)-x(qx,2))*l2(jB)<0);
        qx = qx((p(jB,1)-x(qx,1))*l1(jA)+(p(jB,2)-x(qx,2))*l2(jA)<0);
        
        H(qx) = l2(jA).*l1(jB)-l1(jA).*l2(jB); % (-1)*curvature
    
    end
    
end

H(H>=0) = 1;
H(H< 0) =-1;

end