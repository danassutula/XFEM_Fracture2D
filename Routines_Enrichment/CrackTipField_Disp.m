function [u,ux,uy] = CrackTipField_Disp(K1,K2,G,k,r,theta)

n  = length(r);
u  = zeros(n,2);

ur = zeros(n,2);
ut = zeros(n,2);
ux = zeros(n,2);
uy = zeros(n,2);

drdx = cos(theta);
drdy = sin(theta);

dtdx = -drdy./r;
dtdy =  drdx./r;

theta = theta/2;
c = cos(theta);
s = sin(theta);
r = sqrt(r);

if K1 ~= 0
    
    const = (K1/2/G/sqrt(2*pi));
    
    u(:,1) = c.*(k-1+s.^2*2).*r*const;
    u(:,2) = s.*(k+1-c.^2*2).*r*const;
    
    const = const*0.5;
    
    ur(:,1) = c.*(k-1+s.^2*2)./r*const;
    ur(:,2) = s.*(k+1-c.^2*2)./r*const;
    
    ut(:,1) = s.*(c.^2*6-k-1).*r*const;
    ut(:,2) = c.*(k+5-c.^2*6).*r*const;
    
end

if K2 ~= 0
    
    const = (K2/2/G/sqrt(2*pi));
    
    u(:,1) = u(:,1) + s.*(k+1+c.^2*2).*r*const;
    u(:,2) = u(:,2) - c.*(k-1-s.^2*2).*r*const;
    
    const = const*0.5;
    
    ur(:,1) = ur(:,1) + s.*(k+1+c.^2*2)./r*const;
    ur(:,2) = ur(:,2) + c.*(1-k+s.^2*2)./r*const;
    
    ut(:,1) = ut(:,1) + c.*(k-3+c.^2*6).*r*const;
    ut(:,2) = ut(:,2) + s.*(k-3+c.^2*6).*r*const;
    
end

for i = 1:n
   
    ux(i,:) = ur(i,:)*drdx(i) + ut(i,:)*dtdx(i);
    uy(i,:) = ur(i,:)*drdy(i) + ut(i,:)*dtdy(i);
    
end

end