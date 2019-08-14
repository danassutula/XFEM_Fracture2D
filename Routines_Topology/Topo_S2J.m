function [mNCrdJ,mLNodJ] = Topo_S2J(mNCrdS,mLNodS,jLNodJ)

mLNodJ = mLNodS(:,jLNodJ);
jNdJ2S = unique(mLNodJ); % already sorted
mNCrdJ = mNCrdS(jNdJ2S,:);

for iNdJac = 1:length(jNdJ2S)
    mLNodJ(mLNodJ == jNdJ2S(iNdJac)) = iNdJac;
end

end