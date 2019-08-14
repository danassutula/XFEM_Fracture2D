function [d1S,d2S] = GrdCrackLen(idtip,istip,x_crk,dx_elm,x_elm,x_crs,q_egx)

%==========================================================================

% for element nodal variations in space
% and for crack tip variations in space

%==========================================================================

d1S = 0;
d2S = 0;

if size(x_crs,1) ~= 2
    warning('no crack-element interaction --> RETURN');  return
end

if all(q_egx == 0)
    warning('no crack-element intersection --> RETURN'); return
end

% (!) The crack segment should always point in direction of tip

if idtip == 1 % so reverse
    x_crk = x_crk([2,1],:);
    x_crs = x_crs([2,1],:);
    q_egx = q_egx([2,1],1);
end

if q_egx(1) ~= 0 % sgm. inwards - subtract variations
    
    [d1s,d2s] = GetRates_edg(istip,x_crs(1,:),x_crk,q_egx(1),dx_elm,x_elm);
    
    d1S = -d1s;
    d2S = -d2s;
    
end

if q_egx(2) ~= 0 % sgm. outwards - add variations (REPETITION)
    
    [d1s,d2s] = GetRates_edg(istip,x_crs(2,:),x_crk,q_egx(2),dx_elm,x_elm);
    
    d1S = d1S + d1s;
    d2S = d2S + d2s;
    
end
end

function [ds,d2s] = GetRates_edg(wCkRot,x_crs,x_crk,i_edg,dx_elm,x_elm)

%==========================================================================

% (!) NOTE: d2x_egX ~= [-dx_egX(2);dx_egX(1)],
%     unless the crack rotates with the edge
%     use: dx = Nj*dxj; d2x = Nj*d2xj + dNj*dxj

% % if crack is non-static / init.
% d2x_egX = [-dx_egX(2);dx_egX(1)];

%==========================================================================

X=Gauss_Glb2Lcl(x_crs,x_elm); [N,D]=ShapesStd_omg(1,X);
dx_egX = (N*dx_elm)'; d2x_egX = [-dx_egX(2);dx_egX(1)];

if wCkRot
    dx_ckX = [x_crk(3)-x_crs(2);x_crs(1)-x_crk(1)];
else
    dx_ckX = [0;0];
end

if i_edg < size(x_elm,1)
    p = [i_edg,i_edg+1];
else
    p = [i_edg,1];
end

x_edg  =  x_elm(p,:);
dx_edg = dx_elm(p,:);

lc = x_crk(2,:)-x_crk(1,:);
le = x_edg(2,:)-x_edg(1,:);

ut = sqrt(lc*lc(:))\lc; dun =-ut; 
un = [ -ut(2), ut(1)];  dut = un; 

if ~wCkRot % no rotation
    dun(:) = 0; dut(:) = 0; 
end

dle = dx_edg(2,:)-dx_edg(1,:);
dle = dle(:); le = le(:); 

xt = ut*le; dxt = dut*le + ut*dle;
xn = un*le; dxn = dun*le + un*dle;

cot_b = xt/xn; % (!?) du == 0
dct_b = (dxt*xn-xt*dxn)/xn^2;

% 1st ORDER VARIATIONS:

% due to element rotation
ds = ut*dx_egX-(un*dx_egX)*cot_b;
% due to crack rotation
ds = ds + (un*dx_ckX)*cot_b;

if ~wCkRot % (!!!) d2x = Nj*d2xj + (dx_crk-dx_edg)*iJ*dNjdX*dxj
    d2x_egX = d2x_egX + ((ut*ds-dx_egX')*((D*x_elm)\(D*dx_elm)))';
end

% 2st ORDER VARIATIONS:

% due to element rot. (n.b. d(un*dx_egX) ~= 0)
d2s = dut*dx_egX - (dun*dx_egX)*cot_b + ...
ut*d2x_egX-(un*d2x_egX)*cot_b-(un*dx_egX)*dct_b;
% due to crack rot. (n.b. d(un*dx_ckX) = 0)
d2s = d2s + (un*dx_ckX)*dct_b;

% % due to crack rotation (complete)
% d2x_ckX = [-dx_ckX(2);dx_ckX(1)];
% d2s = d2s + (dun*dx_ckX)*cot_b + ...
% (un*d2x_ckX)*cot_b+(un*dx_ckX)*dct_b;

end