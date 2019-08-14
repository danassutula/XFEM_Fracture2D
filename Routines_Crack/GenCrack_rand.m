function [nCrack,cCkCrd] = GenCrack_rand(n_crack,crack_len,crack_sep,box_xlim,box_ylim,FOS_crack2box)
%
% Random Crack Generator

if nargin == 5
    FOS_crack2box = 0;
end

crack_len = crack_len/2; % crack half-length
crack_sep = crack_sep^2; % (easier to compare)

scale = 1e5;

it = round(rand(n_crack,1)*scale); it(it==0) = 1;
ix = round(rand(n_crack,1)*scale); ix(ix==0) = 1;
iy = round(rand(n_crack,1)*scale); iy(iy==0) = 1;

t = linspace(0,pi,scale);
x = linspace(box_xlim(1),box_xlim(2),scale);
y = linspace(box_ylim(1),box_ylim(2),scale);

t = t(it);
x = x(ix);
y = y(iy);

id = true(n_crack,1);

% remove close-by cracks
for i = 1:n_crack-1
    
    xi = x(i);
    yi = y(i);
    
    for j = i+1:n_crack
        
        if (x(j)-xi)^2+(y(j)-yi)^2 < crack_sep
           id(j) = false; 
        end
        
    end
end

if FOS_crack2box > 0
    
    % remove cracks close to box
    
    box_xlim(1) = box_xlim(1)+crack_len*FOS_crack2box;
    box_xlim(2) = box_xlim(2)-crack_len*FOS_crack2box;
    
    box_ylim(1) = box_ylim(1)+crack_len*FOS_crack2box;
    box_ylim(2) = box_ylim(2)-crack_len*FOS_crack2box;
    
    for i = 1:n_crack
        
        if x(i) < box_xlim(1) || x(i) > box_xlim(2) || ...
           y(i) < box_ylim(1) || y(i) > box_ylim(2)
            
            id(i) = false; % remove crack, too close to boundary
            
        end
        
    end
    
end

id = find(id);

t = t(id);
x = x(id);
y = y(id);

nCrack = length(id);
cCkCrd = cell(nCrack,1);

for i = 1:nCrack
    
   n = [cos(t(i)), sin(t(i))];
   T = [-n;n]; xc = [x(i),y(i)];
   
   cCkCrd{i} = xc([1,1],:) + crack_len*T;
   
end

