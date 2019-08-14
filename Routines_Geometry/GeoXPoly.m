function [X,ID] = GeoXPoly(p1,p2,N)

% N = num. of intersection point desired

n1 = length(p1);
n2 = length(p2);

if nargin == 2
    N = n1+n2; % estimate
end

X  = zeros(N,2); 
ID = zeros(N,2);

j0 = [0,1];
k  = 0;

for i1 = 1:n1-1; j1 = j0+i1;
    for i2 = 1:n2-1; j2 = j0+i2;
        
        x = GeoXLine(p1(j1,:),p2(j2,:));
        
        if ~isempty(x); k=k+1; 
            
            X(k,1) = x(1);
            X(k,2) = x(2);
            
            ID(k,1) = j1(1); 
            ID(k,2) = j2(1);
            
            if k == N
                break
            end
            
        end
        
    end
    
    if k == N
        break
    end
    
end

k = 1:k; X = X(k,:); ID = ID(k,:);

end