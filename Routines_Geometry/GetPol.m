function [t,r] = GetPol(x,y)

%==========================================================================
% Take precaution when coordinate is on one of the crack faces, i.e.:
% x = -ev & y = 0 ==> t = pi or -pi
%==========================================================================

r = sqrt(x*x+y*y); % ~x10 faster than using norm

if r == 0
    t = 0;
else
    t = acos(x/r);
    
    if y < 0
        t = -t;
    end
end

%==========================================================================

end