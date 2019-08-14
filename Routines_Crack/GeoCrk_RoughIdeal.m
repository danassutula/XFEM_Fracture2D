function [Rq,Rv,Rp,Rt,xi,yi,y0] = GeoCrk_RoughIdeal(cCkCrd)
%
% Compute ideal roughnesses of a crack distribution given by 'cCkCrd'.
% The cracks MUST be: x-sequential, non-overlapping and non-intersecting.
% The roughnesses are computed by simply connecting the crack tips.
%
% NOTE! When comparing the results to those obtained by using the function
% 'GeoCrk_Rough' it needs to be born in mind that 'GeoCrk_Rough' computes
% roughness only of cracks that have two intersections; hence, the end
% cracks are ignored. In other words, 'IdealRoughness' and 'GeoCrk_Rough'
% will give different results, unless the cracks used in 'IdealRoughness'
% is 'cCkCrd = cCkCrd(2:end-1)'.

%--------------------------------------------------------------------------
% Get Profile
%--------------------------------------------------------------------------

xi = cell2mat(cCkCrd);

yi = xi(:,2);
xi = xi(:,1);

%--------------------------------------------------------------------------
% Get roughness
%--------------------------------------------------------------------------
if length(xi)>1
    
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
