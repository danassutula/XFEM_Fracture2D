function [b,db] = Growth_KinkAngle_localSym_Williams(r)

% Second order accurate small kink angle approximation.
% Equivalent to the hoop stress criterion \cite{Karihaloo1982}.


% r = SIF ratio, K2/K1
% b = crack tip kink angle.

tol = 1e-6;

b = max(-2*abs(r),-acos(1/3)); % [O(b^2),b_sup]
b = b*sign(r);

db = inf;

while abs(db) > tol

    % C11 =  0.25*(3*cos(0.5*b)+cos(1.5*b));
    % C12 = -0.75*(sin(0.5*b)+sin(1.5*b));
    
    C21 =  0.25*(sin(0.5*b)+sin(1.5*b));
    C22 =  0.25*(cos(0.5*b)+3*cos(1.5*b));
    
    dC21 = 0.125*cos(0.5*b)+0.375*cos(1.5*b);
    dC22 =-0.125*sin(0.5*b)-1.125*sin(1.5*b);
    
    db =-(C21+C22*r)/(dC21+dC22*r);

    b = b + db;
    
end