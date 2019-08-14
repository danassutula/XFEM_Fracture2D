function [sxx,syy,sxy] = CrackTipField_Stress(K1,K2,r,theta)

n = size(r);
sxx = zeros(n);
syy = zeros(n);
sxy = zeros(n);

theta = theta/2;
c1  = cos(theta);
s1  = sin(theta);

theta = theta*3;
c3 = cos(theta);
s3 = sin(theta);

if K1 ~= 0
    const = (K1/sqrt(2*pi))./sqrt(r);
    sxx = const.*c1.*(1-s1.*s3);
    syy = const.*c1.*(1+s1.*s3);
    sxy = const.*c1.*s1.*c3;
end

if K2 ~= 0
    const = (K2/sqrt(2*pi))./sqrt(r);
    sxx = sxx - const.*s1.*(2+c1.*c3);
    syy = syy + const.*s1.*c1.*c3;
    sxy = sxy + const.*c1.*(1-s1.*s3);
end

end