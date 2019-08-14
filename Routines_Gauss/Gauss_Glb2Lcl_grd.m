function X = Gauss_Glb2Lcl_grd(x,xn,dxn)

error('This one is not used. Use Gauss_Glb2Lcl instead')

% warning message: g.p. mapping did not converge
msg = 'mapping Gauss points; residual, dX = ';

nx = size(x,1);
nn = size(xn,1);

if ~(nx && nn)
    X = []; return
end

% control
n_max = 3; tol = 1e-12; tol = tol^2;

% shift coord. to center (better accuracy when mapping)
xm = sum(xn,1)/nn;

xn(:,1) = xn(:,1) - xm(1);
xn(:,2) = xn(:,2) - xm(2);

x(:,1) =  x(:,1) - xm(1);
x(:,2) =  x(:,2) - xm(2);

X(nx,2) = 0; dXdx(2,2) = 0; Xi_zrs(1,2) = 0;

for ix = 1:nx
    
    n = 0; xi = x(ix,:); Xi = Xi_zrs; dSi = 1;
    
    while  dSi>tol && n<n_max
        
        [N,dNdX] = LgBasis_omg(Xi,nn);
        
        dxdX = dNdX*xn;
        detJ = dxdX(1)*dxdX(4)-dxdX(2)*dxdX(3);
        
        dXdx(1) =  dxdX(4);
        dXdx(2) = -dxdX(2);
        dXdx(3) = -dxdX(3);
        dXdx(4) =  dxdX(1);
        
        dXi = detJ\(xi-N*xn)*dXdx;
        dSi = dXi(1)^2+dXi(2)^2;
        
        Xi = Xi+dXi; n = n+1;
        
    end
    
    if n==n_max && dSi>tol
        warning([msg,num2str(sqrt(dSi))])
    end
    
    X(ix,:) = Xi;
    
end
end
