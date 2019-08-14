function bool = GeoElm_xTip(x_tip,x_elm)

% bool == 0, if point is outside
% bool == 1, if point is inside

h = (x_elm(2,1)-x_elm(1,1))^2 ... 
    + (x_elm(2,2)-x_elm(1,2))^2;

tol = 1e-9*sqrt(h);

if x_tip(1)+tol < min(x_elm(:,1)) || x_tip(1)-tol > max(x_elm(:,1)) || ...
   x_tip(2)+tol < min(x_elm(:,2)) || x_tip(2)-tol > max(x_elm(:,2))
    
    bool = 0;
    
else
    
    bool = 1;
    
    n = size(x_elm,1);r = zeros(n,2);
    r(:,1) = x_tip(1,1) - x_elm(:,1);
    r(:,2) = x_tip(1,2) - x_elm(:,2);
    
    tol = -tol^2;
    
    for q = [1:n;2:n,1]
        if  r(q(1),1)*r(q(2),2)-r(q(2),1)*r(q(1),2) < tol
            bool = 0; break
        end
    end
    
end
end