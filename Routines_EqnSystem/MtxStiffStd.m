function Kg = MtxStiffStd(nGlDof,mNdCrd,mElNod,vElPhz,cDMatx,dNdX,w)

nDimes = 2;
nGauss = length(w);

[nElems,nElNod] = size(mElNod);
mElDof = zeros(nDimes,nElNod);

nElDof = nElNod*nDimes;
nElmLK = nElDof^2;
nSprGK = nElmLK*nElems;

if nSprGK == 0
    Kg = sparse([],[],[],nGlDof,nGlDof);
    return
end

Ke0 = zeros(nElDof,nElDof);

jSprLK = 1:nElmLK;
vIdGrd = ones(1,nElDof);

vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

B  = zeros(3,nElDof);
iB = nDimes:nDimes:nElDof;
jB = iB - 1;

igd0 = 1:nDimes;
igp0 = 1:nGauss;
dXdx = zeros(2);

vElPhz = vElPhz(:)';

for iPhase = 1:length(cDMatx)

    % constit. mtx. for phase
    D = cDMatx{iPhase};
    
for iElemn = find(vElPhz==iPhase)
    
    mElCrd = mNdCrd(mElNod(iElemn,:),:);
    
    %{
    vLCrd0 = sum(mElCrd,1)/nElNod;
    mElCrd(:,1) = mElCrd(:,1)-vLCrd0(1);
    mElCrd(:,2) = mElCrd(:,2)-vLCrd0(2);
    %}
    
    %----------------------------------------------------------------------
    % Element Stiffness Matrix
    %----------------------------------------------------------------------
    
    Ke = Ke0;
    igd = igd0;
    
    for igp = igp0;
        
        J = dNdX(igd,:)*mElCrd;
        
        %------------------------------------------------------------------
        % Strain Matrix (2D specific)
        %------------------------------------------------------------------
        detJ = J(1)*J(4)-J(2)*J(3);
        
        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
        
        dNdx = dXdx*dNdX(igd,:);
        
        B(1,jB) = dNdx(1,:);
        B(3,iB) = dNdx(1,:);
        B(2,iB) = dNdx(2,:);
        B(3,jB) = dNdx(2,:);
        %------------------------------------------------------------------
        
        Ke = Ke + B'*D*B.*(detJ*w(igp));
        
        igd = igd + nDimes;
        
    end
    %----------------------------------------------------------------------
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    mElDof(1,:) = nDimes*mElNod(iElemn,:)+(1-nDimes);
    for i = 2:nDimes; mElDof(i,:) = mElDof(i-1,:) + 1; end
    
    vElDof = mElDof(:);
    mRwGrd = vElDof(:,vIdGrd);
    mClGrd = mRwGrd';
    
    jSprRw(jSprLK) = mRwGrd(:);
    jSprCl(jSprLK) = mClGrd(:);
    vSprGK(jSprLK) = Ke(:);
    
    jSprLK = jSprLK + nElmLK;
    %----------------------------------------------------------------------
    
end
end

Kg = sparse(jSprRw,jSprCl,vSprGK,nGlDof,nGlDof);

end
