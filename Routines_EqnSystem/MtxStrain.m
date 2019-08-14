function mBMatx = MtxStrain(nElNod,mCrDrv)

%--------------------------------------------------------------------------
% Standard strain matrix for 2D
% nElNod = number of nodal variables
%--------------------------------------------------------------------------

% nDimes = 2;

mBMatx = zeros(3,2*nElNod);

i = 1;
j = 2;

for iElNod = 1:nElNod
    
    mBMatx(1,i) = mCrDrv(1,iElNod);
    mBMatx(2,j) = mCrDrv(2,iElNod);
    mBMatx(3,i) = mCrDrv(2,iElNod);
    mBMatx(3,j) = mCrDrv(1,iElNod);
    
    i = i + 2;
    j = j + 2;
    
end
end
