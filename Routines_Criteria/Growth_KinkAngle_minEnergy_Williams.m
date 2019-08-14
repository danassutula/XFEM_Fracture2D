function [b,db] = Growth_KinkAngle_minEnergy_Williams(r)

% Second order accurate small kink angle approximation.
% Equivalent to the hoop stress criterion \cite{Karihaloo1982}.

% r = SIF ratio, K2/K1
% b = crack tip kink angle.

tol = 1e-6;

b = max(-2*abs(r),-acos(1/3)); % [O(b^2),b_sup]
b = b*sign(r);

db = inf;

while abs(db) > tol

    C11 =  0.25*(3*cos(0.5*b)+cos(1.5*b));
    C12 = -0.75*(sin(0.5*b)+sin(1.5*b));
    C21 =  0.25*(sin(0.5*b)+sin(1.5*b));
    C22 =  0.25*(cos(0.5*b)+3*cos(1.5*b));
    
    dC11 = -0.375*sin(0.5*b)-0.375*sin(1.5*b);
    dC12 = -0.375*cos(0.5*b)-1.125*cos(1.5*b);
    dC21 =  0.125*cos(0.5*b)+0.375*cos(1.5*b);
    dC22 = -0.125*sin(0.5*b)-1.125*sin(1.5*b);
    
    d2C11 = -0.1875*cos(0.5*b)-0.5625*cos(1.5*b);
    d2C12 =  0.1875*sin(0.5*b)+1.6875*sin(1.5*b);
    d2C21 = -0.0625*sin(0.5*b)-0.5625*sin(1.5*b);
    d2C22 = -0.0625*cos(0.5*b)-1.6875*cos(1.5*b);
    
    % ke  = (C11+C12*r)^2+(C21+C22*r)^2;
    dke  = 2*((C11+C12*r)*(dC11+dC12*r) +(C21+C22*r)*(dC21+dC22*r));
    d2ke = 2*((C11+C12*r)*(d2C11+d2C12*r)+(C21+C22*r)*(d2C21+d2C22*r) + ...
        (dC11+dC12*r)^2 +(dC21+dC22*r)^2);
    
    db  = -dke/d2ke;

    b = b + db;
    
end