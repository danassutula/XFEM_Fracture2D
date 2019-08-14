function [u,ux] = CrackTipField_Jump(K1,K2,G,k,r)

K = (k+1)/G/sqrt(2*pi)*[K2,K1];
r = sqrt(r);

u  = r*K;
ux = -(0.5./r)*K; % x is in opposite direction to r

end