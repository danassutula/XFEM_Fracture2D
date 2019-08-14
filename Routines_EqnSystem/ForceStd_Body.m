function Fg = ForceStd_Body(nGlDof,mNdCrd,mElNod,N,dNdX,w,vElLod)

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

N_zro = zeros(1,nElNod);

for iElemn = 1:nElems
    
    %----------------------------------------------------------------------
    % Element Force Vector
    %----------------------------------------------------------------------

    mElCrd = mNdCrd(mElNod(iElemn,:),:);

    N_sum = N_zro;
    i_gdv = j_gdv;
    
    for i = j_gpt
        
        J = dNdX(i_gdv,:)*mElCrd;
        detJ = J(1)*J(4)-J(2)*J(3);
        
        N_sum = N_sum + N(i,:).*(detJ*w(i));
        i_gdv = i_gdv + nDimes;
        
    end
    
    % 2D
    vLForc(jElDof+1) = N_sum*vElLod(1);
    vLForc(jElDof+2) = N_sum*vElLod(2);
    
    %----------------------------------------------------------------------
    % Global Assembly (Sparse)
    %----------------------------------------------------------------------

    mElDof(1,:) = nDimes*mElNod(iElemn,:)+(1-nDimes);
    for i = 2:nDimes; mElDof(i,:) = mElDof(i-1,:) + 1; end

    jSprRw(jSprLF) = mElDof(:);
    vSprGF(jSprLF) = vLForc;

    jSprLF = jSprLF + nElDof;
    
end

Fg = sparse(jSprRw,ones(nSprGF,1),vSprGF,nGlDof,1);

end
