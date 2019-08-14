function [b_ref,KK_ref] = Growth_DiscreteAngle_energy

% n. discrete eval.'s
n = 1e3+1;

% analytical (theta_max = 75.6 , Gc_I = 1, Gc_II = 2.546)
% numerical  (theta_max = 75.22, Gc_I = 1, Gc_II = 2.546) (okey)
t_max = 75.22*pi/180;

h_ref = t_max/(n-1);
b_ref = linspace(-t_max-h_ref,h_ref,n+2);

% G = ( A11*K1^2 + 2*A12*K1*K2 + A22*K2^2)/E_str
A11 = 4*((1-b_ref/pi)./(1+b_ref/pi)).^(b_ref/pi)./(4-sin(b_ref).^2).^2;
A12 = A11.*(-2*sin(2*b_ref));
A22 = A11.*(4+5*sin(b_ref).^2);
A11 = A11.*(4-3*sin(b_ref).^2);

% numerical diff.
h_ref = 2*h_ref; % x2

a11 = (A11(3:end)-A11(1:end-2))/h_ref;
a12 = (A12(3:end)-A12(1:end-2))/h_ref;
a22 = (A22(3:end)-A22(1:end-2))/h_ref;

% A11 = A11(2:end-1);
% A12 = A12(2:end-1);
% A22 = A22(2:end-1);

% K_II / K_I
KK_ref = (-a12-sqrt(a12.^2-a11.*a22))./a22; KK_ref(end) = 0;
b_ref = b_ref(2:end-1);

end