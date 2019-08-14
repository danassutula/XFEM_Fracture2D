function [mNoInc,cCkCrd,mTpAct,mTpRdi,mCk2Up] = ...
    GeoCrk_Intersect_boundary_frz(mCk2Up,cCkCrd,mTpAct,mTpRdi,f_xrs)

%--------------------------------------------------------------------------
% Check for intersections on domain boundaries (BC); freeze tip position
%--------------------------------------------------------------------------

global mNdCrd 
global cFxNod cLdNod
global cBnNod nBound

mNoInc = zeros(size(mCk2Up));

for i = find(mCk2Up(:,1) & mTpAct(:,1))'
    for j = 1:nBound
        
        [~,s,~,j_sgm] = LevelSet2Poly_bnd(cCkCrd{i}(1,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,1)*f_xrs
            
            bool = 0;
            
            for k = 1:length(cFxNod)
                if any(cFxNod{k} == cBnNod{j}(j_sgm,1))
                    bool = 1; break
                end
            end
            
            if ~bool
                for k = 1:length(cLdNod)
                    if any(cLdNod{k}(:,1) == cBnNod{j}(j_sgm,1))
                        bool = 1; break
                    end
                end
            end
            
            if bool
                
                cCkCrd{i} = cCkCrd{i}(2:end,:);
                
                mCk2Up(i,1) = 0;
                mTpAct(i,1) = 0;
                mNoInc(i,1) =-1;
                
            end
            
            break
            
        end
    end
end

for i = find(mCk2Up(:,2) & mTpAct(:,2))'
    for j = 1:nBound
        
        [~,s,~,j_sgm] = LevelSet2Poly_bnd(cCkCrd{i}(end,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,2)*f_xrs
            
            bool = 0;
            
            for k = 1:length(cFxNod)
                if any(cFxNod{k} == cBnNod{j}(j_sgm,1))
                    bool = 1; break
                end
            end
            
            if ~bool
                for k = 1:length(cLdNod)
                    if any(cLdNod{k}(:,1) == cBnNod{j}(j_sgm,1))
                        bool = 1; break
                    end
                end
            end
            
            if bool
                
                cCkCrd{i} = cCkCrd{i}(1:end-1,:);
                
                mCk2Up(i,2) = 0;
                mTpAct(i,2) = 0;
                mNoInc(i,2) =-1;
                
            end
            
            break
            
        end
    end
end
end
