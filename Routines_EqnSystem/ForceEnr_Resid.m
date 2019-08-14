function Fg = ForceEnr_Resid(nGlDof,mNdCrd,mLNodS,vElEnr,cLNodE,mNDofE,...
                             cGsDvS,cGsDvE,cGsWgt,mElLod)

% RHS = - int B*(D*eps_0 + sig_0) dOmg

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
dXdx = zeros(2,2);

for iElemn = vElEnr(:)'
    vElLod = mElLod(:,iElemn);
    
    if any(vElLod)
        
        %------------------------------------------------------------------
        % Element Force Vector
        %------------------------------------------------------------------
        
        mElCrd = mNdCrd(mLNodS(iElemn,:),:);
        
        dNdX_s = cGsDvS{iElemn};
        dNdX_e = cGsDvE{iElemn};
        w      = cGsWgt{iElemn};
        
        vLNodE = cLNodE{iElemn};
        nLNodE = length(vLNodE);
        
        mLDofE = mNDofE(:,vLNodE);
        nLDofE = nLNodE*nDimes;
        
        vLForc = zeros(nLDofE,1);
        
        B_sum = zeros(nDimes,nLNodE);
        i_gdv = i0_gdv;
        
        for i = 1:length(w)
            
            J = dNdX_s(i_gdv,:)*mElCrd;
            % detJ = J(1)*J(4)-J(2)*J(3);
            
            dXdx(1) =  J(4); dXdx(3) = -J(3);
            dXdx(2) = -J(2); dXdx(4) =  J(1); % /detJ
            
            B_sum = B_sum + (dXdx*dNdX_e(i_gdv,:)).*w(i);
            i_gdv = i_gdv + nDimes;
            
        end
        
        % 2D
        vLForc(1:nDimes:nLDofE) = B_sum(1,:)*vElLod(1)+B_sum(2,:)*vElLod(3);
        vLForc(2:nDimes:nLDofE) = B_sum(2,:)*vElLod(2)+B_sum(1,:)*vElLod(3);
        
        %------------------------------------------------------------------
        % Global Assembly (Sparse)
        %------------------------------------------------------------------
        
        jSprLF = jSprLF(end)+(1:nLDofE);
        jSprRw(jSprLF) = mLDofE(:);
        vSprGF(jSprLF) = vLForc;
        
    end
end

j = jSprLF(end)+1:nSprGF;

jSprRw(j) = [];
vSprGF(j) = [];

% do not forget force is negative, i.e. F = -int B*sig_0 dOmg
Fg = sparse(jSprRw,ones(length(jSprRw),1),-vSprGF,nGlDof,1);

end
