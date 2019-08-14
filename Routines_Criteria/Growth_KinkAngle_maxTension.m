function b = Growth_KinkAngle_maxTension(r) % [v (Poisson's ratio)]

% r = SIF ratio, K2/K1
% b = crack tip kink angle.

%%% Maximum surface energy density cirterion:
% theta_c = argmax[ SIGMA_nn(theta)^2 + (1+v)*TAU_nt(theta)^2 ]
%
%%% (equivalent to the hoop-stress criterion)
% theta_c = argmax[ SIGMA_n(theta) ]
%
% Both critria correspond to a plane of zero shear and maximum tension.

% E_n = (K1*cos(b/2)^3-3*K2*sin(b/2)*cos(b/2)^2)^2 + ...
%   (1+v)*(K1*sin(b/2)*cos(b/2)^2 + K2*cos(b/2)*(1-3*sin(b/2)^2))^2;

b = zeros(size(r));

b(r == 0) = 0;
b(r == inf) = -acos(1/3);
q = ~(r == 0 || r == inf);
b(q) = 2*(atan((1-sqrt(1+8*r(q).^2))./(4*r(q))));

end