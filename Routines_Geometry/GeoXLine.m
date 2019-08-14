function x = GeoXLine(v1,v2)

x = []; tol = 1e-9; 

n1 = v1(2,:)-v1(1,:);
n2 = v2(2,:)-v2(1,:);

a = n1(1)*n2(2)-n1(2)*n2(1);
b = v2(2,:)-v1(1,:);

% check if lines not parallel (alternative)
% if a^2 > tol^2*(n1(1)*n1(1)+n1(2)*n1(2))*(n2(1)*n2(1)+n2(2)*n2(2))

if a*a > 0 % check if lines not parallel (seems to work)
    
    % scaling factors
    % l = [n1(:),n2(:)]\b(:)
    
    % scalling factors (faster)
    l1 = (b(1)*n2(2)-b(2)*n2(1))/a;
    l2 = (b(2)*n1(1)-b(1)*n1(2))/a;
    
    if l1>-tol && l1<1+tol && l2>-tol && l2<1+tol
        x = v1(1,:) + l1*n1;
    end
    
end