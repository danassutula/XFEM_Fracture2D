function Fg = ForceEnr_Body(nGlDof,mNdCrd,mLNodS,vElEnr,cLNodE,mNDofE,...
                            cGsDvS,cGsShE,cGsWgt,vElLod)

nDimes = 2;

nElEnr = length(vElEnr);
nLNodS = size(mLNodS,2);

nLDofS = nLNodS*nDimes;
nSprGF = nElEnr*nLDofS*5; % predict size

if nSprGF == 0
    Fg = sparse([],[],[],nGlDof,1);
    return
end

jSprLF = 0;

jSprRw = zeros(nSprGF,1);
vSprGF = zeros(nSprGF,1);

i0_gdv = 1:nDimes;

for iElemn = vElEnr(:)'
    
    %----------------------------------------------------------------------
    % Element Force Vector
    %----------------------------------------------------------------------
    
    mElCrd = mNdCrd(mLNodS(iElemn,:),:);
    
    dNdX = cGsDvS{iElemn};
    N    = cGsShE{iElemn};
    w    = cGsWgt{iElemn};
    
    vLNodE = cLNodE{iElemn};
    nLNodE = length(vLNodE);
    
    mLDofE = mNDofE(:,vLNodE);
    nLDofE = nLNodE*nDimes;
    
    vLForc = zeros(nLDofE,1);
    
    N_sum = zeros(1,nLNodE);
    i_gdv = i0_gdv;
    
    for i = 1:length(w)
        
        J = dNdX(i_gdv,:)*mElCrd;
        detJ = J(1)*J(4)-J(2)*J(3);
        
        N_sum = N_sum + N(i,:).*(detJ*w(i));
        i_gdv = i_gdv + nDimes;
        
    end
    
    % 2D
    vLForc(1:nDimes:nLDofE) = N_sum*vElLod(1);
    vLForc(2:nDimes:nLDofE) = N_sum*vElLod(2);
    
    %------------------------------------------------------------------
    % Global Assembly (Sparse)
    %------------------------------------------------------------------
    
    jSprLF = jSprLF(end)+(1:nLDofE);
    jSprRw(jSprLF) = mLDofE(:);
    vSprGF(jSprLF) = vLForc;
    
end

j = jSprLF(end)+1:nSprGF;

jSprRw(j) = [];
vSprGF(j) = [];

Fg = sparse(jSprRw,ones(length(jSprRw),1),vSprGF,nGlDof,1);

end
