function [f,Df] = EnrJmp_Brn(x,x_tip) % x_elm

% Jump in the Branch function: sin(lambda*theta)

x(:,1) = x(:,1) - x_tip(1);
x(:,2) = x(:,2) - x_tip(2);

r = sqrt(x(:,1).^2+x(:,2).^2);

% 1/2 singularity:
% l = 0.5; jump = 2;

f  =  2.*sqrt(r);
Df = -1./sqrt(r);

%{
% General singularity:
l = 0.45; jump = 2*sin(l*pi);

f  =  jump*r.^l;
Df = -jump*l*r.^(l-1);
%}

end