function [Rq,Rv,Rp,Rt,xi,yi,y0] = GeoCrk_RoughPost(cCkCrd,mCkJun)
%
% Roughness is computed for cracks that have a pair of intersections.
% Therefore, END-CRACKS ARE IGNORED (since only a single intersection).
% If a crack has more than two intersections, the fracture profile becomes
% non-unique (or at least unclear); this results in a warning message.

msg = 'not clear how to compute roughness; crack intersected more than twice';

%--------------------------------------------------------------------------
% Pre-process cracks (i.e. remove blending segments)
%--------------------------------------------------------------------------

nCrack = length(cCkCrd);

for i = 1:nCrack;
    
    if mCkJun(i,1)
        i_bgn = mCkJun(i,1)+1;
    else
        i_bgn = 1;
    end
    
    if mCkJun(i,2)
        i_end = mCkJun(i,2);
    else
        i_end = length(cCkCrd{i});
    end
    
    % no blending segments
    cCkCrd{i} = cCkCrd{i}(i_bgn:i_end,:);
    
end

%--------------------------------------------------------------------------
% Find Crack Intersections
%--------------------------------------------------------------------------

XS = zeros(10*nCrack,2); k_XS = 0;
ID = zeros(nCrack,2);    k_ID = zeros(nCrack,1);

for i1 = 1:nCrack-1
    X1 = cCkCrd{i1};
    
    for i2 = i1+1:nCrack
        X2 = cCkCrd{i2};
        
        [xs,id] = GeoXPoly(X1,X2,1);
        
        if ~isempty(xs)
        
            k_XS     = k_XS     + 1 ; XS(k_XS,:)      = xs;
            k_ID(i1) = k_ID(i1) + 1 ; ID(i1,k_ID(i1)) = id(1);
            k_ID(i2) = k_ID(i2) + 1 ; ID(i2,k_ID(i2)) = id(2);
            
        end
        
    end
end

jCrack = find(k_ID > 1);
nCrack = length(jCrack);

if any(k_ID(jCrack)>2)
   warning(msg)
end

%--------------------------------------------------------------------------
% Get Profile
%--------------------------------------------------------------------------

if nCrack > 0
    
    ID = sort(ID(jCrack,1:2),2);
    XS = XS(1:k_XS,:);
    
    for i = 1:nCrack
        cCkCrd{i} = cCkCrd{jCrack(i)}(ID(i,1)+1:ID(i,2),:);
    end
    
    cCkCrd = [cat(1,cCkCrd{1:nCrack}) ; XS];
    [~,j] = unique(round(cCkCrd(:,1)*1e12));
    
    xi = cCkCrd(j,1);
    yi = cCkCrd(j,2);
    
else
    
    xi  = [];
    yi  = [];
    
end

%--------------------------------------------------------------------------
% Get roughness
%--------------------------------------------------------------------------
if ~isempty(xi)
    
    % profile length
    L = xi(end)-xi(1);
    
    % segments' length
    dx = xi(2:end)-xi(1:end-1);
    
    % neautral axis, i.e. \int y dx = 0
    y0 = sum((yi(2:end)+yi(1:end-1)).*dx)/(2*L);
    
    % shifted profile
    yi = yi - y0;
    
    % RMS roughness
    Rq = sqrt(sum((yi(2:end).^2+yi(2:end).*yi(1:end-1)+yi(1:end-1).^2).*dx)/(3*L));
    
    % maximum valley depth
    Rv = min(yi);
    
    % maximum peak height
    Rp = max(yi);
    
    % maximum height of the profile
    Rt = Rp-Rv;
    
    % restore profile
    yi = yi + y0;
    
else
    Rq=[]; Rv=[]; Rp=[]; Rt=[]; y0=[];
end
