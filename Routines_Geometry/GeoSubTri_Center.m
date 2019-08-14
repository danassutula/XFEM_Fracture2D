function lnod = GeoSubTri_Center(x)

npt = size(x,1); % subjective, but works good as an estimate
pts = 2:npt;

lnod = zeros(10*npt,3);
lnod(:,1) = 1;

nel = 0; % for counting delaunay elements

x(:,1) = x(:,1) - x(1,1); 
x(:,2) = x(:,2) - x(1,2); 

a = zeros(2,2);
b = zeros(2,2);

for i = 2:(npt-1)
    a(1,:) = x(i,:);
    
    for j = (i+1):npt
        a(2,:) = x(j,:);
        
        bool = 1;
        
        for k = pts(~ismember(pts,[i,j]))
            b(2,:) = 2*x(k,:);
            
            if ~isempty(GeoXLine(a,b)) % ture of el. crosses another el.
                bool = 0;
                break
            end
            
        end
        
        if bool
            nel = nel + 1;
            lnod(nel,2) = i;
            lnod(nel,3) = j;
        end
        
    end
end

lnod = lnod(1:nel,:);

for i = 1:nel
    
    v1 = x(lnod(i,2),:) - x(lnod(i,1),:);
    v2 = x(lnod(i,3),:) - x(lnod(i,2),:);
    
    if v1(1)*v2(2)-v1(2)*v2(1) < 0
        lnod(i,:) = lnod(i,[1,3,2]);
    end
    
end
