function nrmLS = GeoCrk_PathComp(cCkCrd_1,mCkJun_1,cCkCrd_2,mCkJun_2)
%
% Compare two fracture paths by computing the L2-norm of the minimum 
% absolute distance between them. The profile given by cCkCrd_1 is used
% as a reference for computing the integral.

nCrack = length(cCkCrd_1);
if nCrack ~= length(cCkCrd_1);
    error('Expected the same number of cracks')
end

if isempty(mCkJun_1)
    mCkJun_1 = zeros(nCrack,2);
end
if isempty(mCkJun_2)
    mCkJun_2 = zeros(nCrack,2);
end
    
nGauss = 5;
[X,W] = Gauss_DomLin(nGauss);
N = ShapesStd_lin(nGauss,X);

if ~isempty(mCkJun_1) % trim off crack blending segments
    [L_tot,~,cCkCrd_1] = GeoCrk_Length(cCkCrd_1,mCkJun_1); else
    [L_tot,~,cCkCrd_1] = GeoCrk_Length(cCkCrd_1,zeros(nCrack,2));
end    

if ~isempty(mCkJun_2) % trim off crack blending segments
    [~,~,cCkCrd_2] = GeoCrk_Length(cCkCrd_2,mCkJun_2); else
    [~,~,cCkCrd_2] = GeoCrk_Length(cCkCrd_2,zeros(nCrack,2));
end   

nrmLS = 0; % init. L2 norm of the level-set between two cracks

for i = 1:nCrack
    for j = 1:size(cCkCrd_1{i},1)-1
        
        L_ths = norm(cCkCrd_1{i}(j+1,:)-cCkCrd_1{i}(j,:));
        
        for k = 1:nGauss
            
            absLS = LevelSet2Poly_abs(N(k,:)*cCkCrd_1{i}([j,j+1],:),cCkCrd_2{i});
            nrmLS = nrmLS + absLS^2*(0.5*L_ths)*W(k);
            
        end
        
    end
end
    
nrmLS = sqrt(nrmLS)/L_tot;
