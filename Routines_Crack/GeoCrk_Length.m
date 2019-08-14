function [l_tot,l_crk,cCkCrd] = GeoCrk_Length(cCkCrd,mCkJun)
%
% l_tot     = combined crack length
% l_crk     = individual crack lengths
% cCkCrd    = trimmed crack coordiantes (i.e. disregarding blending branches) 

n = length(cCkCrd);
l_crk = zeros(n,1);

for i = 1:n
    
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
    
    % combined crack length
    l_crk(i) = sum(sqrt(...
        (cCkCrd{i}(2:end,1)-cCkCrd{i}(1:end-1,1)).^2 + ...
        (cCkCrd{i}(2:end,2)-cCkCrd{i}(1:end-1,2)).^2));
    
end

% combined length of all cracks
l_tot = sum(l_crk);
