function y_bar = GeoCrk_MeanDepth(cCkCrd)
%
% get mean depth of crack distribution
% each crack is assumed to be straight

nCrack = length(cCkCrd);

y = 0;
l = 0;

for i = 1:nCrack
    
    d = norm(diff(cCkCrd{i}([1,end],1)));
    y = y + mean(cCkCrd{i}([1,end],2))*d;
    l = l + d;
    
end

y_bar = y/l;
