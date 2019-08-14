function [mFnVal,mFnDrv,mNdShf] = EnrFun_Brn(mElCrd,mGsCrd,nEnFun,tpCrd,tpAlf)
%==========================================================================

% NOTE with regard to the enrichment function:
%   atan2(   0,-1) = pi
%   atan2(-eps,-1) =-pi

nElNod = size(mElCrd,1);
nGauss = size(mGsCrd,1);

mFnVal = zeros(nGauss,nEnFun);
mFnDrv = zeros(nGauss*2,nEnFun);
mNdShf = zeros(nEnFun,nElNod);

mGsCrd(:,1) = mGsCrd(:,1)-tpCrd(1);
mGsCrd(:,2) = mGsCrd(:,2)-tpCrd(2);

mElCrd(:,1) = mElCrd(:,1)-tpCrd(1); 
mElCrd(:,2) = mElCrd(:,2)-tpCrd(2);

s = sin(tpAlf);
c = cos(tpAlf);
T = [c,-s;s,c];

mGsCrd = mGsCrd*T; % transform to local axis (crack dir.)
mElCrd = mElCrd*T; % ...

[theta,r] = cart2pol(mGsCrd(:,1),mGsCrd(:,2));
[mFnVal,mFnDrv] = Branch(r,theta); % here, mFnDrv are local

j = 1:2;

% transoform local to global
for i = 1:nGauss
    mFnDrv(j,:) = T*mFnDrv(j,:); j = j + 2;
end

[theta,r] = cart2pol(mElCrd(:,1),mElCrd(:,2));
mNdShf = Branch(r,theta)';

%==========================================================================

end

function [f,Df] = Branch(r,theta)

%==========================================================================
% DETAIL
%{
F1 = sqrt(r)*sin(t/2);
F2 = sqrt(r)*cos(t/2);
F3 = sqrt(r)*sin(t/2)*cos(t/2)^2;
F4 = sqrt(r)*cos(t/2)*sin(t/2)^2;

F1x' = -0.5*sin(t/2)/sqrt(r);
F2x' =  0.5*cos(t/2)/sqrt(r);
F3x' = -0.5*sin(t/2)*cos(t/2)^2*(4*cos(t/2)^2-3)/sqrt(r);
F4x' = -0.5*cos(t/2)*sin(t/2)^2*(4*cos(t/2)^2-1)/sqrt(r);

F1y' = 0.5*cos(t/2)/sqrt(r);
F2y' = 0.5*sin(t/2)/sqrt(r);
F3y' = 0.5*cos(t/2)*(4*cos(t/2)^4-5*cos(t/2)^2+2)/sqrt(r);
F4y' = 0.5*sin(t/2)*(4*cos(t/2)^4-3*cos(t/2)^2+1)/sqrt(r);

,where x' & y' are in local cord. sys.

[dx'/dx, dy'/dx; dx'/dy, dy'/dy] = [cos(a), -sin(a); sin(a), cos(a)]
%}
%==========================================================================

% make data in a column order
r = r(:); theta = theta(:);

%==========================================================================
% Sqrt singularity, l = 1/2
%
%   enrichment functions for Ux & Uy:
%   \phi = sqrt(r)*[sin(t/2),cos(t/2),sin(t/2)*sin(t),cos(t/2)*sin(t)]
%==========================================================================

n = 2*length(r);

s = sin(theta/2);
c = cos(theta/2);

S = 2*s.*c;
C = c.^2-s.^2;

Sc = S.*c;
Ss = S.*s;

sqrt_r = sqrt(r);

f = repmat(sqrt_r,1,4).*[s,c,Ss,Sc];

sqrt_r = repmat(0.5./sqrt_r,1,4);

Df = zeros(n,4);
Df(1:2:n-1,:) = sqrt_r.*[C.*s-Sc,C.*c+Ss,-S.*Sc-C.*Ss,S.*Ss-C.*Sc];
Df(2:2:n  ,:) = sqrt_r.*[C.*c+Ss,Sc-C.*s,(1+C.*C).*s+C.*Sc,(1+C.*C).*c-C.*Ss];

% f  =  f(:,1:3);
% Df = Df(:,1:3);

%==========================================================================
% General singularity, l
%
%   enrichment functions for Ux & Uy:
%   \phi = r^l*[sin(l*t),cos(l*t),...
%               sin(l*t)*sin(t)^2,cos(l*t)*cos(t)^2,...
%               sin(l*t)*sin(2*t),cos(l*t)*sin(2*t)]
%==========================================================================

%{
l = 0.5;

J(2,2)  = 0;
f(1,6)  = 0;
Df(2,6) = 0;

S = sin(l*theta);
C = cos(l*theta);

s = sin(theta);
c = cos(theta);

f(1) = S;
f(2) = C;
f(3) = S*s^2;
f(4) = C*c^2;
f(5) = S*s*c;
f(6) = C*s*c;

J(1) = c;
J(2) = s;
J(3) =-s/r;
J(4) = c/r;

Df(1,:) = (r^(l-1)*l).*f;

f = r^l.*f;

Df(2,1) = l*C;
Df(2,2) =-l*S;
Df(2,3) = l*C*s^2 + 2*S*s*c;
Df(2,4) =-l*S*c^2 - 2*C*s*c;
Df(2,5) = l*C*s*c + S*c^2 - S*s^2;
Df(2,6) =-l*S*s*c + C*c^2 - C*s^2;

Df(2,:) = r^l.*Df(2,:);

Df = J*Df;

% f = f(1:5);
% Df = Df(:,1:5);
%}
     
end
