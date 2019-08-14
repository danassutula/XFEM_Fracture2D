function Gs = IntJ(cCkCrd,mNDspS,mNDspE,P0,phz,D,R)

%--------------------------------------------------------------------------
% Crack Tip Energy Release Rate: The J-Integral
%   1. No pressure,
%--------------------------------------------------------------------------

% Using a transformation matrix:
% T = [  c^2,  s^2,  2*s*c; ...
%        s^2,  c^2, -2*s*c; ...
%       -s*c,  s*c,  c^2-s^2 ];

global mNdCrd
global mLNodS
global cLNodE

global mGsStd_omgDrv
global vGsStd_omgWgt

global cGsEnr_omgDvS
global cGsEnr_omgDvE
global cGsEnr_omgDvJ
global cGsEnr_omgWgt

nDimes = 2;

% energy release rate
Gs = zeros(length(cCkCrd),2);

% transformations
T2 = zeros(2,2);
T3 = zeros(3,3);

dXdx = zeros(2,2);
e    = zeros(3,1);
i0_gdv = 1:nDimes;

for i_crk = 1:length(cCkCrd)
    
    mTpCrd = cCkCrd{i_crk}([1,end],:);
    
    % dir. needs to always point towards crack tip
    mTpDir(1,:) = cCkCrd{i_crk}(1,:)   - cCkCrd{i_crk}(2,:);
    mTpDir(2,:) = cCkCrd{i_crk}(end,:) - cCkCrd{i_crk}(end-1,:);
    
    for i_tip = 1:2
        if R(i_crk,i_tip) > 0
            
            Gs_tip = 0;
            
            vTpCrd = mTpCrd(i_tip,:);
            vTpDir = mTpDir(i_tip,:);
            
            n = vTpDir/sqrt(vTpDir(1)^2+vTpDir(2)^2);
            
            T2(1) = n(1); T2(3) =-n(2);
            T2(2) = n(2); T2(4) = n(1);
            
            T3(1) =    n(1)^2;  T3(4) =    n(2)^2;  T3(7) =   2*n(1)*n(2); 
            T3(2) =    n(2)^2;  T3(5) =    n(1)^2;  T3(8) =  -2*n(1)*n(2);
            T3(3) =-n(1)*n(2);  T3(6) = n(1)*n(2);  T3(9) = n(2)^2-n(1)^2;
            
            % detect elements within domain; get smoothing weight
            [vDmElm,mDmWgt] = Elems_sif(mNdCrd,mLNodS,vTpCrd,R(i_crk,i_tip));
            
            for i_dom = 1:length(vDmElm)
                i_elm = vDmElm(i_dom);
                
                % align with tip
                mLDspS = mNDspS(:,mLNodS(i_elm,:))'*T2;
                
                % good to set relative to tip
                mElCrd = mNdCrd(mLNodS(i_elm,:),:);
                mElCrd(:,1) = mElCrd(:,1)-vTpCrd(1);
                mElCrd(:,2) = mElCrd(:,2)-vTpCrd(2);
                
                % align element with tip
                mElCrd = mElCrd*T2;
                
                % rotate pre-load and get pre-strain
                D_elm = D{phz(i_elm)};
                e0 = D_elm\(T3*P0(:,i_elm));
                
                % weight function
                vDmWgt = mDmWgt(:,i_dom);
                
                if ~isempty(cLNodE{i_elm})
                    
                    % align with tip
                    mLDspE = mNDspE(:,cLNodE{i_elm})'*T2;
                    
                    dudX = cGsEnr_omgDvS{i_elm}*mLDspS + ...
                           cGsEnr_omgDvE{i_elm}*mLDspE;
                    dxdX = cGsEnr_omgDvS{i_elm}*mElCrd;
                    dqdX = cGsEnr_omgDvS{i_elm}*vDmWgt;
                    
                    w = cGsEnr_omgWgt{i_elm};
                    
                else
                    
                    dudX = mGsStd_omgDrv*mLDspS;
                    dxdX = mGsStd_omgDrv*mElCrd;
                    dqdX = mGsStd_omgDrv*vDmWgt;
                    
                    w = vGsStd_omgWgt;
                    
                end
                
                i_gdv = i0_gdv;
                
                for i = 1:length(w);
                    
                    J = dxdX(i_gdv,:); detJ = det(J);
                    
                    % already aligned with tip - no transformation.
                    dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                    dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                    
                    % already aligned with tip - no transformation.
                    ux = dXdx*dudX(i_gdv,:);
                    qx = dXdx*dqdX(i_gdv);
                    
                    % get strain
                    e(1) = ux(1);
                    e(2) = ux(4);
                    e(3) = ux(2)+ux(3);
                    
                    % add pre-strain
                    e = e + e0;
                    
                    % get total stress
                    s = D_elm*e;
                    
                    % external work
                    W = (ux(1,1)*s(1)+ux(1,2)*s(3))*qx(1) + ...
                        (ux(1,1)*s(3)+ux(1,2)*s(2))*qx(2);
                    
                    % internal energy
                    U = 0.5*(e(1)*s(1)+e(2)*s(2)+e(3)*s(3))*qx(1);
                    
                    % energy release rate
                    Gs_tip = Gs_tip + (W - U)*detJ*w(i);
                    
                    i_gdv = i_gdv + nDimes;
                    
                end
                
            end
            
            Gs(i_crk,i_tip) = Gs_tip;
            
        end
    end
end
end