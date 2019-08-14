function [mNoInc,cCkCrd,mTpAct,mTpRdi,mCkJun] = ...
    GeoCrk_Intersect_boundary(mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,f_xrs)

%--------------------------------------------------------------------------
% Check for crack intersections with domain boundaries
%--------------------------------------------------------------------------

global mNdCrd cBnNod;
nBound=length(cBnNod);

mNoInc = zeros(size(mCk2Up));

for i = find(mCk2Up(:,1) & mTpAct(:,1))'
    for j = 1:nBound
        
        [~,s,d] = LevelSet2Poly_bnd(cCkCrd{i}(1,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,1)*f_xrs
            if (cCkCrd{i}(1,:)-cCkCrd{i}(2,:))*d(:) > 0
                
                cCkCrd{i} = [(1+1e-4)*d+cCkCrd{i}(1,:);cCkCrd{i}];
                
                mTpRdi(i,1) = 0;
                mTpAct(i,1) = 0;
                mNoInc(i,1) = 1;
                
                if mCkJun(i,2)
                    mCkJun(i,2) = mCkJun(i,2)+1;
                end
                
                break
                
            end
        end
    end
end

for i = find(mCk2Up(:,2) & mTpAct(:,2))'
    for j = 1:nBound
        
        [~,s,d] = LevelSet2Poly_bnd(cCkCrd{i}(end,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,2)*f_xrs
            if (cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:))*d(:) > 0
                
                cCkCrd{i} = [cCkCrd{i};cCkCrd{i}(end,:)+(1+1e-4)*d];
                
                mTpRdi(i,2) = 0;
                mTpAct(i,2) = 0;
                mNoInc(i,2) = 1;
                
                break
                
            end
        end
    end
end
