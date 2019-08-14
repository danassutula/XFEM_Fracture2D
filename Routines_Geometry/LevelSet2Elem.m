function [x_lvl,s_lvl] = LevelSet2Elem(x_pnt,x_elm)

%--------------------------------------------------------------------------

n = size(x_elm,1);
q = [1:n;2:n,1];

s_lvl = inf;
x_lvl = [ ];

for i = 1:n

    x_lin = x_elm(q(:,i),:);
    
    t1 = x_lin(2,1)-x_lin(1,1);
    t2 = x_lin(2,2)-x_lin(1,2);
    
    r1 = x_pnt(1)-x_lin(1,1);
    r2 = x_pnt(2)-x_lin(1,2);
    
    l = sqrt(t1*t1+t2*t2);
    t1 = t1/l; t2 = t2/l;
    s = r1*t1+r2*t2;
    
    if s < 0
        x = x_lin(1,:);
    elseif s > l
        x = x_lin(2,:);
    else
        x = x_lin(1,:)+[s*t1,s*t2];
    end
    
    s = (x(1)-x_pnt(1))^2 + (x(2)-x_pnt(2))^2;
    
    if s_lvl > s
        s_lvl = s;
        x_lvl = x;
    end

end

s_lvl = sqrt(s_lvl);

end