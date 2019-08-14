function [mNdCrd,mElNod] = GeoSubTri_Refine(mNdCrd,mElNod)

% only for triangles, i.e. size(mElNod,2) = 3

%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================

nPoint = size(mNdCrd,1);
nElems = size(mElNod,1);

nPtNew = nPoint+nElems;
nElNew = nElems*3;

mNdCrd(nPoint+1:nPtNew,:) = zeros(nElems,2);
mElNod(nElems+1:nElNew,:) = zeros(nElNew-nElems,3);

nl = nElems;
nx = nPoint;

for il = 1:nElems

    vElNod = mElNod(il,:);
    
    nx = nx+1;
    xn = mNdCrd(vElNod,:);
    xc = 3\(xn(1,:)+xn(2,:)+xn(3,:));
    
    mNdCrd(nx,:) = xc;
    
    ln = GeoSubTri_Center([xc;xn]);
    
    mElNod(il,1) = nx;
    mElNod(il,2) = vElNod(ln(1,2)-1);
    mElNod(il,3) = vElNod(ln(1,3)-1);
    
    nl = nl+1;
    mElNod(nl,1) = nx;
    mElNod(nl,2) = vElNod(ln(2,2)-1);
    mElNod(nl,3) = vElNod(ln(2,3)-1);
    
    nl = nl+1;
    mElNod(nl,1) = nx;
    mElNod(nl,2) = vElNod(ln(3,2)-1);
    mElNod(nl,3) = vElNod(ln(3,3)-1);
    
end

%==========================================================================

end