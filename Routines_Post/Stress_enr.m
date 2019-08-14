function [s_ij,s_xy] = Stress_enr(mNDspS,mNDspE,mNdCrd,vElEnr,mLNodS,cLNodE,...
                                  vElPhz,nPhase,cDMatx,mPrLod,...
                                  cGsShS,cGsDvS,cGsDvE)

%--------------------------------------------------------------------------
% Stress at Gauss points (enr. el.)
%--------------------------------------------------------------------------

vElEnr = vElEnr(:)';

ngt = 10e6; % estimate

s_ij = zeros(3,ngt);
s_xy = zeros(ngt,2);
e    = zeros(3,1);

dXdx = zeros(2,2);
igd0 = 1:2;
igt  = 0;

for i_phz = 1:nPhase
    
    % constit. mtx. for phase
    D = cDMatx{i_phz};
    
for iel = vElEnr(vElPhz(vElEnr)==i_phz)
    
    mElCrd = mNdCrd(mLNodS(iel,:),:);
    
    %{
    vLCrd0 = sum(mElCrd,1)/nLNodS;
    mElCrd(:,1) = mElCrd(:,1)-vLCrd0(1);
    mElCrd(:,2) = mElCrd(:,2)-vLCrd0(2);
    %}
    
    mLDspS = mNDspS(:,mLNodS(iel,:))';
    mLDspE = mNDspE(:,cLNodE{iel})';
    
    x    = cGsShS{iel}*mElCrd;
    dxdX = cGsDvS{iel}*mElCrd;
    dudX = cGsDvS{iel}*mLDspS+cGsDvE{iel}*mLDspE;
    
    ngp = size(x,1);
    
    s_xy(igt+1:igt+ngp,:) = x;
    s0 = mPrLod(:,iel);
    
    igd = igd0;
    
    for igp = 1:ngp
        
        J = dxdX(igd,:); detJ = J(1)*J(4)-J(2)*J(3);
        
        dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
        dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
        
        dudx = dXdx*dudX(igd,:);
        
        igd = igd + 2;
        igt = igt + 1;
        
        e(1) = dudx(1);
        e(2) = dudx(4);
        e(3) = dudx(2)+dudx(3);
        
        s_ij(:,igt) = D*e + s0;
        
    end
    
end
end

s_xy = s_xy(1:igt,:);
s_ij = s_ij(:,1:igt);

end