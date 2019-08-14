function Fg = ForceStd_Resid(nGlDof,mNdCrd,mElNod,dNdX,w,mElLod)

% RHS = - int B*(D*eps_0 + sig_0) dOmg

nDimes = 2;
nGauss = length(w);

[nElems,nElNod] = size(mElNod);
mElDof = zeros(nDimes,nElNod);

nElDof = nElNod*nDimes;
vLForc = zeros(nElDof,1);

nSprGF = nElDof*nElems;
jSprRw = zeros(nSprGF,1);
vSprGF = zeros(nSprGF,1);

jSprLF = 1:nElDof;
jElDof = 0:nDimes:(nElDof-nDimes);

j_gpt = 1:nGauss;
j_gdv = 1:nDimes;
dXdx = zeros(2,2);

B_zro = zeros(nDimes,nElNod);

for iElemn = 1:nElems
    vElLod = mElLod(:,iElemn);
    
    if any(vElLod)
        
        %------------------------------------------------------------------
        % Element Force Vector
        %------------------------------------------------------------------
        
        mElCrd = mNdCrd(mElNod(iElemn,:),:);
        
        B_sum = B_zro;
        i_gdv = j_gdv;
        
        for i = j_gpt
            
            J = dNdX(i_gdv,:)*mElCrd;
            % detJ = J(1)*J(4)-J(2)*J(3);
            
            dXdx(1) =  J(4); dXdx(3) = -J(3);
            dXdx(2) = -J(2); dXdx(4) =  J(1); % /detJ
            
            B_sum = B_sum + (dXdx*dNdX(i_gdv,:)).*w(i); % *detJ
            i_gdv = i_gdv + nDimes;
            
        end
        
        % 2D
        vLForc(jElDof+1) = B_sum(1,:)*vElLod(1)+B_sum(2,:)*vElLod(3);
        vLForc(jElDof+2) = B_sum(2,:)*vElLod(2)+B_sum(1,:)*vElLod(3);
        
    else
        
        vLForc(:) = 0;
        
    end
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------
    
    mElDof(1,:) = nDimes*mElNod(iElemn,:)+(1-nDimes);
    for i = 2:nDimes; mElDof(i,:) = mElDof(i-1,:) + 1; end
    
    jSprRw(jSprLF) = mElDof(:);
    vSprGF(jSprLF) = vLForc;
    
    jSprLF = jSprLF + nElDof;
    
    %----------------------------------------------------------------------
    
end

% do not forget force is negative, i.e. F = -int B*sig_0 dOmg
Fg = sparse(jSprRw,ones(nSprGF,1),-vSprGF,nGlDof,1);

end
