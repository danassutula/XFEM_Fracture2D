function Gs_mix = IntM(nCrack,cCkCrd,mNDspS,mNDspE,P0,phz,D,G,k,R,K1,K2)

% Can include:
%   Gs_mix
%   Gs_st1
%   Gs_st2
%   Gs_com

%==========================================================================
% Interaction Integral (M): element contour form of integration
%   1. No pressure, 
%   2. No body force,
%   3. No thermal force.
%==========================================================================

%--------------------------------------------------------------------------
% Define Global Variables
%--------------------------------------------------------------------------
global mNdCrd
global mLNodS
global cLNodE

global mGsStd_omgShp
global mGsStd_omgDrv
global vGsStd_omgWgt

global cGsEnr_omgShS
global cGsEnr_omgDvS
global cGsEnr_omgDvE
global cGsEnr_omgWgt

global cBrStd_eTp
%--------------------------------------------------------------------------

%==========================================================================
%                               BULLETPROOF!!!
%==========================================================================

nDimes = 2;

Gs_mix = zeros(nCrack,2);
% Gs_st1 = zeros(nCrack,2);
% Gs_st2 = zeros(nCrack,2);
% Gs_com = zeros(nCrack,2);

% std. state
s = zeros(2,2);
% aux. state
S = zeros(2,2);

% transformations
T2 = zeros(2,2);
T3 = zeros(3,3);

dXdx = zeros(2,2);

i0_gdv = 1:nDimes;

for i_crk = 1:nCrack
    
    mTpCrd = cCkCrd{i_crk}([1,end],:);
    
    % dir. needs to always point towards crack tip
    mTpDir(1,:) = cCkCrd{i_crk}(1,:)   - cCkCrd{i_crk}(2,:);
    mTpDir(2,:) = cCkCrd{i_crk}(end,:) - cCkCrd{i_crk}(end-1,:);
    
    for i_tip = 1:2
        if R(i_crk,i_tip) > 0
            
            gs_mix = 0;
%             gs_st1 = 0;
%             gs_st2 = 0;
            
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
            
            % only used in the analytical solution (assume constant)
            i = phz(cBrStd_eTp{i_crk,i_tip}(1));
            G_tip = G(i); k_tip = k(i);
            
            for i_dom = 1:length(vDmElm)
                i_elm = vDmElm(i_dom);
                
                % align with tip
                mLDspS = mNDspS(:,mLNodS(i_elm,:))'*T2;
                
                % good to set relative to tip
                mElCrd = mNdCrd(mLNodS(i_elm,:),:);
                mElCrd(:,1) = mElCrd(:,1)-vTpCrd(1);
                mElCrd(:,2) = mElCrd(:,2)-vTpCrd(2);
                
                % align with tip
                mElCrd = mElCrd*T2;
                
                % get Hook for el.
                D_elm = D{phz(i_elm)};
                C_elm = inv(D_elm);
                
                % rotate pre-load and get pre-strain
                e0 = C_elm*(T3*P0(:,i_elm));
                
                % weight function
                vDmWgt = mDmWgt(:,i_dom);
                
                if ~isempty(cLNodE{i_elm})
                    
                    % align with tip
                    mLDspE = mNDspE(:,cLNodE{i_elm})'*T2;
                    
                    dudX = cGsEnr_omgDvS{i_elm}*mLDspS + ...
                           cGsEnr_omgDvE{i_elm}*mLDspE;
                    dxdX = cGsEnr_omgDvS{i_elm}*mElCrd;
                    dqdX = cGsEnr_omgDvS{i_elm}*vDmWgt;
                    
                    x = cGsEnr_omgShS{i_elm}*mElCrd;
                    w = cGsEnr_omgWgt{i_elm};
                    
                else
                    
                    dudX = mGsStd_omgDrv*mLDspS;
                    dxdX = mGsStd_omgDrv*mElCrd;
                    dqdX = mGsStd_omgDrv*vDmWgt;
                    
                    x = mGsStd_omgShp*mElCrd;
                    w = vGsStd_omgWgt;
                    
                end
                
                %----------------------------------------------------------
                % State 2 (auxiliary)
                %----------------------------------------------------------
                [t,r] = cart2pol(x(:,1),x(:,2));
                
                [~,DUDX] = CrackTipField_Disp(K1,K2,G_tip,k_tip,r,t);
                [SXX,SYY,SXY] = CrackTipField_Stress(K1,K2,r,t);
                
                EXX = C_elm(1)*SXX + C_elm(4)*SYY + C_elm(7)*SXY;
                EYY = C_elm(2)*SXX + C_elm(5)*SYY + C_elm(8)*SXY;
                EXY = C_elm(3)*SXX + C_elm(6)*SYY + C_elm(9)*SXY;
                %----------------------------------------------------------
                
                i_gdv = i0_gdv;
                
                for i = 1:length(w)
                    
                    J = dxdX(i_gdv,:); detJ = det(J); dA = detJ*w(i);
                    
                    % already aligned with tip. No transformation.
                    dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                    dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                    
                    % already aligned with tip. No transformation.
                    ux = dXdx*dudX(i_gdv,:);
                    qx = dXdx*dqdX(i_gdv);
                    
                    exx =       ux(1) + e0(1);
                    eyy =       ux(4) + e0(2);
                    exy = ux(2)+ux(3) + e0(3);
                    
                    sxx = D_elm(1)*exx + D_elm(4)*eyy + D_elm(7)*exy;
                    syy = D_elm(2)*exx + D_elm(5)*eyy + D_elm(8)*exy;
                    sxy = D_elm(3)*exx + D_elm(6)*eyy + D_elm(9)*exy;
                    
                    % Get Cauchy (std)
                    s(1) = sxx; s(3) = sxy;
                    s(2) = sxy; s(4) = syy;
                    
                    % Get Cauchy (aux)
                    S(1) = SXX(i); S(3) = SXY(i);
                    S(2) = SXY(i); S(4) = SYY(i);
                    
                    %------------------------------------------------------
                    %  Mixed State
                    %------------------------------------------------------
                    W = (ux(1,:)*S+DUDX(i,:)*s)*qx;
                    U = 0.5*(exx*SXX(i) + eyy*SYY(i) + exy*SXY(i) + ...
                             EXX(i)*sxx + EYY(i)*syy + EXY(i)*sxy)*qx(1);
                    
                    gs_mix = gs_mix + (W - U)*dA;
                    %------------------------------------------------------   
                    
%                     %------------------------------------------------------
%                     % State 1
%                     %------------------------------------------------------
%                     W = ux(1,:)*s*qx;
%                     U = 0.5*(exx*sxx+eyy*syy+exy*sxy)*qx(1);
%                     
%                     gs_st1 = gs_st1 + (W - U)*dA;
%                     %------------------------------------------------------
                    
%                     %------------------------------------------------------
%                     % State 2 (auxiliary)
%                     %------------------------------------------------------
%                     W = DUDX(i,:)*S*qx;
%                     U = 0.5*(EXX(i)*SXX(i)+EYY(i)*SYY(i)+EXY(i)*SXY(i))*qx(1);
% 
%                     gs_st2 = gs_st2 + (W - U)*dA;
%                     %------------------------------------------------------
                    
                    i_gdv = i_gdv + nDimes;
                    
                end
                
                
            end
            
            Gs_mix(i_crk,i_tip) = gs_mix;
%             Gs_st1(i_crk,i_tip) = gs_st1;
%             Gs_st2(i_crk,i_tip) = gs_st2;
            
        end
    end
end

%==========================================================================

end