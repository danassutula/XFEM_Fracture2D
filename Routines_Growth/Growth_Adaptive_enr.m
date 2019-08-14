
%==========================================================================
% Adaptive crack tip enrichment radius
%==========================================================================


%--------------------------------------------------------------------------
if 0 % based on the crack tip element size
%--------------------------------------------------------------------------

for i = find(mCk2Up(:,1) & mTpAct(:,1))'
    if ~isempty(cBrStd_eTp{i,1}) 
        
        x_elm  = mNdCrd(mLNodS(cBrStd_eTp{i,1}(1),:),:);
        he = ElemSize(x_elm); mTpRdi(i,1) = he*f_tip;
        
    end
end

for i = find(mCk2Up(:,2) & mTpAct(:,2))'
    if ~isempty(cBrStd_eTp{i,2})
        
        x_elm  = mNdCrd(mLNodS(cBrStd_eTp{i,2}(1),:),:);
        he = ElemSize(x_elm); mTpRdi(i,2) = he*f_tip;
        
    end
end

%--------------------------------------------------------------------------
else % based on the crack tip median element size
%-------------------------------------------------------------------------- 

for i = find(mCk2Up(:,1) & mTpAct(:,1))'
    if ~isempty(cBrBln_eFl{i,1})
        
        q   = cBrBln_eFl{i,1};    % crack tip el. id
        tmp = zeros(length(q),1); % to hold el. size
        
        for j = 1:length(q)
            tmp(j) = ElemSize(mNdCrd(mLNodS(q(j),:),:));
        end
        
        mTpRdi(i,1) = f_tip*median(tmp);

    end
end

for i = find(mCk2Up(:,2) & mTpAct(:,2))'
    if ~isempty(cBrBln_eFl{i,2})
        
        q   = cBrBln_eFl{i,2};    % crack tip el. id
        tmp = zeros(length(q),1); % to hold el. size
        
        for j = 1:length(q)
            tmp(j) = ElemSize(mNdCrd(mLNodS(q(j),:),:));
        end
        
        mTpRdi(i,2) = f_tip*median(tmp);

    end
end
end
