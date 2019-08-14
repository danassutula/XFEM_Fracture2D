function [k1,k2] = Growth_GetKinkedSIF_op2(K1,K2,b)

% Hayashi and Nemat-Nasser, 1981
% Energy-release rate and crack kinking under combined loading

% BACKGROUND:
% Based on the maximum energy-release-rate criterion, kinking from a
% straight crack is investigated under the plane strain condition.
% Solutions are obtained by a method that models a kink as a continuous
% distribution of edge dislocations. The energy release rate is expressed
% as a quadratic from of the stress-intensity factors that exist prior to
% the onset of kinking, and the coefficients of this quadratic form are
% tabulated for various values of the kink angle. The examination of the
% results shows that Irwin's formula for the energy releas rate remains
% valid for any kink angle provided that the stress intensity factors in
% the formula are taken equal to those existing at the tip of a vanishingly
% small kink.

% KI  = KI1  * k1 + KI2  * k2
% KII = KII1 * k1 + KII2 * k2
% w   = theta/pi

w   = [0.0 -0.04    -0.08    -0.12    -0.16    -0.20    -0.24    -0.28    -0.32    -0.36    -0.40    -0.44    -0.48    -0.52    -0.56    -0.60    -0.64    -0.68    -0.72    -0.76    -0.80   ];
K11 = [1.0  0.99410  0.97655  0.94794  0.90913  0.86127  0.80579  0.74427  0.67837  0.60981  0.54024  0.47126  0.40426  0.34049  0.28095  0.22632  0.17664  0.12248  0.10537  0.07269  0.04716];
K12 = [0.0  0.18770  0.37069  0.54441  0.70469  0.84784  0.97083  1.07134  1.14784  1.19960  1.22672  1.22996  1.21082  1.17129  1.11390  1.04149  0.95746  0.87185  0.75787  0.65584  0.55040];
K21 = [0.0 -0.06251 -0.12320 -0.18029 -0.23219 -0.27751 -0.31514 -0.34430 -0.36449 -0.37560 -0.37782 -0.37159 -0.35768 -0.33705 -0.31080 -0.28024 -0.24697 -0.21803 -0.17069 -0.13673 -0.10410];
K22 = [1.0  0.98772  0.95131  0.89211  0.81224  0.71460  0.60262  0.48016  0.35137  0.22048  0.09165 -0.03116 -0.14440 -0.24495 -0.33023 -0.39818 -0.44724 -0.47259 -0.49145 -0.48255 -0.45761];

k1 = zeros(size(b));
k2 = k1;

w = pi*w;

for i = 1:numel(b)
    
    if b(i) < 0
        h = 1;
    else
        h =-1;
    end
    
    b_i = b(i)*h;
    
    k1(i) = interp1(w,K11.*K1(i)+h*K12.*K2(i),b_i,'spline');
    k2(i) = interp1(w,h*K21.*K1(i)+K22.*K2(i),b_i,'spline');
    
end

end