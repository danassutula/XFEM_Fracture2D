function Kg = MtxStiffEnr(nGlDof,mNdCrd,mLNodS,vElEnr,cLNodE,mNDofE,...
                          vElPhz,nPhase,cDMatx,cGsDvS,cGsDvE,cGsWgt)

nDimes = 2;

vElEnr = vElEnr(:)';
nElEnr = length(vElEnr);

nLNodS = size(mLNodS,2);
mLDofS = zeros(nDimes,nLNodS);

nLDofS = nLNodS*nDimes;
nSprGK = nElEnr*(nLDofS*5)^2; % predict size

if nSprGK == 0
    Kg = sparse([],[],[],nGlDof,nGlDof);
    return
end

vIdGdS = ones(1,nLDofS);
jSprLK = 0;

vSprGK = zeros(nSprGK,1);
jSprRw = zeros(nSprGK,1);
jSprCl = zeros(nSprGK,1);

Bs  = zeros(3,nLDofS);
iBs = nDimes:nDimes:nLDofS;
jBs = iBs - 1;

dXdx = zeros(2);
igd0 = 1:nDimes;

for iPhase = 1:nPhase

    % constit. mtx. for phase
    D = cDMatx{iPhase};
    
for uElEnr = vElEnr(vElPhz(vElEnr)==iPhase)
    
    mElCrd = mNdCrd(mLNodS(uElEnr,:),:);
    
    %{
    vLCrd0 = sum(mElCrd,1)/nElNod;
    mElCrd(:,1) = mElCrd(:,1)-vLCrd0(1);
    mElCrd(:,2) = mElCrd(:,2)-vLCrd0(2);
    %}
    
    dNdX_s = cGsDvS{uElEnr};
    dNdX_e = cGsDvE{uElEnr};
    
    w = cGsWgt{uElEnr};
    
    %----------------------------------------------------------------------
    % Element Stiffness Matrix (enriched)
    %----------------------------------------------------------------------
    vLNodE = cLNodE{uElEnr};
    nLNodE = length(vLNodE);
    nLDofE = nLNodE*nDimes;
    
    Be  = zeros(3,nLDofE);
    iBe = nDimes:nDimes:nLDofE;
    jBe = iBe - 1;
    
    Kse = zeros(nLDofS,nLDofE);
    Kee = zeros(nLDofE);
    
    igd = igd0;
    for igp = 1:length(w)
        
        J = dNdX_s(igd,:)*mElCrd;
        
        %------------------------------------------------------------------
        % Strain Matrix (2D specific)
        %------------------------------------------------------------------
        detJ = J(1)*J(4)-J(2)*J(3);
        
        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
        
        dNdx_s = dXdx*dNdX_s(igd,:);
        dNdx_e = dXdx*dNdX_e(igd,:);
        %------------------------------------------------------------------
        % B_std
        %------------------------------------------------------------------
        Bs(1,jBs) = dNdx_s(1,:);
        Bs(3,iBs) = dNdx_s(1,:);
        Bs(2,iBs) = dNdx_s(2,:);
        Bs(3,jBs) = dNdx_s(2,:);
        %------------------------------------------------------------------
        % B_enr
        %------------------------------------------------------------------
        Be(1,jBe) = dNdx_e(1,:);
        Be(3,iBe) = dNdx_e(1,:);
        Be(2,iBe) = dNdx_e(2,:);
        Be(3,jBe) = dNdx_e(2,:);
        %------------------------------------------------------------------
        
        DBJw = D*Be.*(detJ*w(igp));
        
        Kse = Kse + Bs'*DBJw;
        Kee = Kee + Be'*DBJw;
        
        igd = igd + nDimes;
        
    end
    %----------------------------------------------------------------------
        
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    mLDofS(1,:) = nDimes*mLNodS(uElEnr,:)+(1-nDimes);
    for i = 2:nDimes; mLDofS(i,:) = mLDofS(i-1,:) + 1; end
    
    mLDofE = mNDofE(:,vLNodE);
    
    vLDofS = mLDofS(:);
    vLDofE = mLDofE(:);
    
    nElmLK = nLDofE*(2*nLDofS+nLDofE);
    jSprLK = jSprLK(end)+(1:nElmLK);
    vIdGdE = ones(1,nLDofE);
    
    mRowSE = vLDofS(:,vIdGdE); mColSE = vLDofE(:,vIdGdS)';
    mRowEE = vLDofE(:,vIdGdE); mColEE = mRowEE';
    
    jSprRw(jSprLK) = [mRowSE(:);mColSE(:);mRowEE(:)];
    jSprCl(jSprLK) = [mColSE(:);mRowSE(:);mColEE(:)];
    vSprGK(jSprLK) = [Kse(:);Kse(:);Kee(:)];
    %----------------------------------------------------------------------
    
end
end

j = jSprLK(end)+1:nSprGK;

jSprRw(j) = [];
jSprCl(j) = [];
vSprGK(j) = [];
    
Kg = sparse(jSprRw,jSprCl,vSprGK,nGlDof,nGlDof);

end
