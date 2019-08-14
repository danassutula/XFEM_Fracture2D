function [b_ref,KK_ref] = Growth_DiscreteAngle_stress

% n. discrete eval.'s
n = 1e3+1;

t_max = 2*acos(sqrt(2/3)); % analytical solution (=70.53deg.)
b_ref = linspace(-t_max,0,n);

% common
x = cos(b_ref/2);

% K_II / K_I
KK_ref = x.*sqrt(1-x.^2)./(x.^2*3-2);

end