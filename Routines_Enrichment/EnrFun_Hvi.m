function [H_gsp,H_nod] = EnrFun_Hvi(mElCrd,mGsCrd,mHvCrd)
%==========================================================================

global tol_ref; % to determine if element node lies on crack

if size(mHvCrd,1) == 2
    
    H_gsp = LevelSet2Line(mGsCrd,mHvCrd); % (signed-distance)
    
    H_gsp(H_gsp >=0) = 1;
    H_gsp(H_gsp < 0) =-1;
    
    H_nod(1,:) = LevelSet2Line(mElCrd,mHvCrd); % (signed-distance)
    H_nod(abs(H_nod)<tol_ref) = 1; % if node lies on line segment
    
    H_nod(H_nod > 0) = 1;
    H_nod(H_nod < 0) =-1;
    
else % nXsCrd == 3 (crack kink)
    
    H_gsp = LevelSet2Poly(mGsCrd,mHvCrd); % (sign)
    [H_nod(1,:),s] = LevelSet2Poly(mElCrd,mHvCrd);
    
    % if node is on poly-line
    H_nod(s<tol_ref) = 1;
    
end
