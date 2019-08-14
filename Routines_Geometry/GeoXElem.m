function [N,X,I] = GeoXElem(x_lin,x_elm)

[N,X,I] = GeoElm_xEdg(x_lin,x_elm);

if N == 1
    
    tol = ((x_elm(2,1)-x_elm(1,1))^2 ...
        +(x_elm(2,2)-x_elm(1,2))^2)*1e-18;

    if GeoElm_xTip(x_lin(1,:),x_elm)
    
        X = [x_lin(1,:);X];
        I = [ 0; I]; N = 2;
        
        if (X(2,1)-X(1,1))^2+(X(2,2)-X(1,2))^2 < tol
            X = X(1,:); I = I(2); N = 1; % redundant
        end
            
    elseif GeoElm_xTip(x_lin(2,:),x_elm)
        
        X = [X;x_lin(2,:)];
        I = [ I; 0]; N = 2;
        
        if (X(2,1)-X(1,1))^2+(X(2,2)-X(1,2))^2 < tol
            X = X(2,:); I = I(1); N = 1; % redundant
        end
        
    end
    
elseif N == 0
    
    if GeoElm_xTip(x_lin(1,:),x_elm)
        X = x_lin; I = [0;0]; N = 2; 
    end
    
end