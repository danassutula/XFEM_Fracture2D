function [N,X,I] = GeoElm_xEdg(x_lin,x_elm)

tol = ((x_elm(2,1)-x_elm(1,1))^2 ...
    + (x_elm(2,2)-x_elm(1,2))^2)*1e-18;

N = 0;
X = zeros(2,2);
I = zeros(2,1);

n = size(x_elm,1);
for i = [1:n;2:n,1]
    x = GeoXLine(x_lin,x_elm(i,:));
    if ~isempty(x); N = N + 1;
        
        X(N,:) = x;
        I(N) = i(1);
        
        if N == 2
            if (X(2,1)-X(1,1))^2 + (X(2,2)-X(1,2))^2 > tol % line on node
                break
            else
                N = N - 1;
            end
        end
    end
end

if N == 2
    if (x_lin(1,1)-X(1,1))^2 + (x_lin(1,2)-X(1,2))^2 > ...
       (x_lin(1,1)-X(2,1))^2 + (x_lin(1,2)-X(2,2))^2
        
        % in succession
        X = X([2,1],:);
        I = I([2,1],:);
    
    end
else
    X = X(1:N,:); % empty if N == 0
    I = I(1:N);
end