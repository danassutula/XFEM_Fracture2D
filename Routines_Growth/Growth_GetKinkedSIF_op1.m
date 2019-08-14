function [k1,k2] = Growth_GetKinkedSIF_op1(K1,K2,b)

% Nuismer, 1975
% Energy-release rate and crack kinking under combined loading

c12 = cos(b/2);
c11 = cos(b);
s11 = sin(b);

k1 = 0.5*c12.*(K1.*(1+c11)-3*K2.*s11);
k2 = 0.5*c12.*(K1.*s11+K2.*(3*c11-1));

% K_eqv = sqrt(k1^2+k2^2);
% K_eqv = sqrt(0.125*(c+1).*( ...
%     (K1.*(1+c)-3*K2.*s).^2+(K1.*s+K2.*(3*c-1)).^2));

end