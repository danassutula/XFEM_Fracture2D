function [p,isAnticlockwise] = Topo_sortNod(p)
%
%--------------------------------------------------------------------------
% 
% Topo_sortNod:
%   Ensures that the sorted boundary points follow an anti-clockwise order
%
% INPUT/OUTPUT:
%   p = [x1,y1;x2,y2;...;x1,y1], closed boundary coordinates that must be
%   sorted either in a clockwise or an anti-clockwise order. The output
%   will be always sorted in an anti-clockwise order.
%
%--------------------------------------------------------------------------

isAnticlockwise = 1;

% v = diff(p,1);              % boundary edge vectors
% b = atan2(v(:,2),v(:,1));   % angle rel. to x-axis

b=atan2(p(2:end,2)-p(1:end-1,2),p(2:end,1)-p(1:end-1,1));      
b(b<0)=b(b<0)+2*pi; b=[b;b(1)]; % make b>0 and close loop

d = 0; % sum diff. in angle

for i = 1:length(b)-1
    if b(i) < pi % && b(i) > 0
        
        if b(i+1)-b(i) < pi
           d=d+b(i+1)-b(i);
        else
           d=d+b(i+1)-b(i)-2*pi;
        end
        
    else % b(i) > pi && b(i) < 2*pi
        
        if b(i)-b(i+1) < pi
           d=d+b(i+1)-b(i);
        else
           d=d+b(i+1)-b(i)+2*pi;
        end
        
    end 
end

if abs(abs(d) - 2*pi) > length(b) * 1e-13 
    error('angle sum does not equal 2*pi')
end

if d < 0 % reverse order
    p = p(end:-1:1,:);
    isAnticlockwise=0;
end
end