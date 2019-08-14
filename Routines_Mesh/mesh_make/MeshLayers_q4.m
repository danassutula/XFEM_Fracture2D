function [coord,lnods,bound,elphz] = MeshLayers_q4(xlim,Y,N,h0)

% -------------------------------------------------------------------------
%
% INPUT:
%   xlim    - limits of rectangular domain on the x-axis
%   Y       - vector of layer interfaces; n. layers = length(Y)-1
%   N       - vector of layer refinements (times layer is finer)
%   h0      - reference element size (N acts on h0 to get h_finer)
%
% OUTPUT:
%   coord   - nodal coordinates
%   lnods   - element topology
%   bound   - boundary nodes (edges: bottom, right, top, left)
%   elphz   - enumerated layers
%
% -------------------------------------------------------------------------

% % from coarse to fine
% dH = 2*hyz*(1-1/2^N);
% 
% % from fine to coarse
% dH = 2*hyf*(2^N-1);

coord = []; % (stacking on top)
lnods = [];
elphz = [];

L = xlim(2)-xlim(1);
H = Y(end)-Y(1);

nx0 = round(L/h0);
if mod(nx0,2)==1;
    nx0 = nx0+1;
end

% DO FIRST LAYER

% mesh size in x
nx = nx0*2^N(1);
hx = L/nx;

% layer depth
dY = Y(2)-Y(1);

if N(2) > N(1) % with step-down
    
    % blending mesh size factor
    f = 2*(1-1/2^(N(2)-N(1))); % (coarse to fine)
    
    % blending mesh depth
    dy = dY/round(dY/hx)*f;
    
    % layer depth
    dY = dY-dy;
    
    % mesh layer
    [coord,lnods] = MeshRect_q4(xlim,[Y(1),Y(2)-dy],nx+1,round(dY/hx)+1);
    
    % mesh layer blending region
    [c,l] = StepDown_q4(xlim,[Y(2)-dy,Y(2)],nx+1,N(2)-N(1)+1);
    
    c = c(nx+2:end,:);
    l = l+size(coord,1)-nx-1;
    
    coord = [coord;c];
    lnods = [lnods;l];
    
else % mesh until end of this layer
    
    % mesh layer (mesh blending region, if any, later i.e. line 73)
    [coord,lnods] = MeshRect_q4(xlim,[Y(1),Y(2)],nx+1,round(dY/hx)+1);
    
end

elphz = ones(size(lnods,1),1);

% DO OTHER LAYERS (EXCEPT LAST)

for i = 2:length(N)-1 % for ith layer
    
    % layer depth
    dY = Y(i+1)-Y(i);
    
    if N(i-1) > N(i) % with step-up
        
        % blending mesh size factor
        f = 2*(2^(N(i-1)-N(i))-1); % (fine to coarse)
        
        % blending mesh depth
        dy = dY/round(dY/hx)*f;
        
        % mesh layer blending region (this is "later")
        [c,l] = StepUp_q4(xlim,[Y(i),Y(i)+dy],nx+1,N(i-1)-N(i)+1);
        
        c = c(nx+2:end,:);
        l = l+size(coord,1)-nx-1;
        
        coord = [coord;c];
        lnods = [lnods;l];
        
        elphz = [elphz;i(ones(size(l,1),1))];
        
        % modify elevation for later
        Y(i) = Y(i) + dy;
        
        % modify layer depth also
        dY = dY-dy;
        
    end
    
    nx = nx0*2^N(i);
    hx = L/nx;
    
    if N(i+1) > N(i) % with step-down
        
        % blending mesh size factor
        f = 2*(1-1/2^(N(i+1)-N(i))); % (coarse to fine)
        
        % blending mesh depth
        dy = dY/round(dY/hx)*f;
        
        % layer depth
        dY = dY-dy;
        
        % mesh layer
        [c,l] = MeshRect_q4(xlim,[Y(i),Y(i+1)-dy],nx+1,round(dY/hx)+1);

        c = c(nx+2:end,:);
        l = l+size(coord,1)-nx-1;
        
        coord = [coord;c];
        lnods = [lnods;l];
        
        elphz = [elphz;i(ones(size(l,1),1))];
        
        % mesh layer blending region
        [c,l] = StepDown_q4(xlim,[Y(i+1)-dy,Y(i+1)],nx+1,N(i+1)-N(i)+1);
        
    else % mesh until end of this layer
        
        % mesh layer
        [c,l] = MeshRect_q4(xlim,[Y(i),Y(i+1)],nx+1,round(dY/hx)+1);
        
    end
    
    c = c(nx+2:end,:);
    l = l+size(coord,1)-nx-1;
    
    coord = [coord;c];
    lnods = [lnods;l];
    
    elphz = [elphz;i(ones(size(l,1),1))];
    
end

% DO LAST LAYERS

i = length(N);

% layer depth
dY = Y(i+1)-Y(i);

if N(i-1) > N(i) % step-up mesh size
    
    % blending mesh size factor
    f = 2*(2^(N(i-1)-N(i))-1); % (fine to coarse)
    
    % blending mesh depth
    dy = dY/round(dY/hx)*f;
    
    % mesh layer blending region 
    [c,l] = StepUp_q4(xlim,[Y(i),Y(i)+dy],nx+1,N(i-1)-N(i)+1);
    
    c = c(nx+2:end,:);
    l = l+size(coord,1)-nx-1;
    
    coord = [coord;c];
    lnods = [lnods;l];
    
    elphz = [elphz;i(ones(size(l,1),1))];
    
    % modify elevation
    Y(i) = Y(i) + dy;
    
    % layer depth
    dY = dY-dy;
    
end

nx = nx0*2^N(i);
hx = L/nx;

% mesh layer
[c,l] = MeshRect_q4(xlim,[Y(i),Y(i+1)],nx+1,round(dY/hx)+1);

c = c(nx+2:end,:);
l = l+size(coord,1)-nx-1;

coord = [coord;c];
lnods = [lnods;l];

elphz = [elphz;i(ones(size(l,1),1))];

% GET DOMAIN BOUNDARY NODES

bound = cell(4,1);
tol=max(L,H)*1e-12;

bound{1} = find( abs(coord(:,2)-min(coord(:,2))) < tol );
bound{2} = find( abs(coord(:,1)-max(coord(:,1))) < tol );
bound{3} = find( abs(coord(:,2)-max(coord(:,2))) < tol );
bound{4} = find( abs(coord(:,1)-min(coord(:,1))) < tol );

end