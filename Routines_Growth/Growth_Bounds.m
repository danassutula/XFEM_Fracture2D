
%--------------------------------------------------------------------------
% Check if cracks are within simulation bounds; stop if true
%--------------------------------------------------------------------------

bool = 0;

for i = find(mCk2Up(:,1))'
    if cCkCrd{i}(1,1) < mCkBox(1,1) || cCkCrd{i}(1,1) > mCkBox(2,1) || ...
       cCkCrd{i}(1,2) < mCkBox(1,2) || cCkCrd{i}(1,2) > mCkBox(2,2)
   
        warning(['Crack #',num2str(i),' tip #',num2str(1),...
            ' outside simulation limits']); bool = 1; break
        
    end
end

if bool
   break 
end

for i = find(mCk2Up(:,2))'
    if cCkCrd{i}(end,1) < mCkBox(1,1) || cCkCrd{i}(end,1) > mCkBox(2,1) || ...
       cCkCrd{i}(end,2) < mCkBox(1,2) || cCkCrd{i}(end,2) > mCkBox(2,2)
   
        warning(['Crack #',num2str(i),' tip #',num2str(2),...
            ' outside simulation limits']); bool = 1; break
        
    end
end

if bool
   break 
end

%--------------------------------------------------------------------------