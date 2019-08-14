function Ms = IntM_CrkLod(nCrack,cCkCrd,ngp,N,w,phz,G,k,R,P,K1,K2)

%==========================================================================

% Crack Tip Energy Release Rate: The M-Integral

%==========================================================================

%--------------------------------------------------------------------------
% Define Global Variables
%--------------------------------------------------------------------------

global mNdCrd
global mLNodS

global cBrStd_eTp
global cHvBln_elm
global cHvStd_elm
global cHvStd_sgm

global with_MapTyp
global tol_abs

%--------------------------------------------------------------------------

nLNodS = size(mLNodS,2);

Ms = zeros(nCrack,2);

for i_crk = 1:nCrack
    
    p = P(i_crk,[2,1])'; % [pressure,shear] -> [shear,pressure] (convenient)
    
    if p(1)^2+p(2)^2 > 0 && ~isempty(cHvStd_elm{i_crk})
        
        nCkCrd = size(cCkCrd{i_crk},1);
        jSgCrd = [2,1;nCkCrd-1,nCkCrd];
        vTpSgm = [1;nCkCrd-1]; % (aux.)
        
        for i_tip = 1:2
            if R(i_crk,i_tip) > 0 && ~isempty(cBrStd_eTp{i_crk,i_tip})
                
                % init. interaction int.
                Ms_tip = 0;
                
                % get tip segment
                mSgCrd = cCkCrd{i_crk}(jSgCrd(i_tip,:),:);
                vTpCrd = mSgCrd(2,:);
                
                % get enr. elements that are on the tip segment
                vElEnr = [cBrStd_eTp{i_crk,i_tip};...
                          cHvBln_elm{i_crk,i_tip};...
                          cHvStd_elm{i_crk}(cHvStd_sgm{i_crk}(:,i_tip) == vTpSgm(i_tip))];
                
                % detect domain elements and get smoothing weight
                [jElDom,jElBnd,mBnWgt] = Elems_tip(...
                    mNdCrd,mLNodS(vElEnr,:),vTpCrd,R(i_crk,i_tip));
                
                % good to transpose here
                mBnWgt = mBnWgt';
                
                % only used in the analytical solution (assume constant)
                i = phz(cBrStd_eTp{i_crk,i_tip}(1));
                G_tip = G(i); k_tip = k(i);
                
                %----------------------------------------------------------
                % Loop over inner-domain elements
                %----------------------------------------------------------
                
                for i_dom = jElDom(:)'
                    
                    mElCrd = mNdCrd(mLNodS(vElEnr(i_dom),:),:);
                    
                    % get intersection of crack with el.
                    [nXsCrd,mXsCrd] = GeoXElem(mSgCrd,mElCrd);
                    
                    % to integrate, must have a length and must have a jump 
                    if nXsCrd == 2 && ~all(LevelSet2Line(mElCrd,mSgCrd)>-tol_abs)
                        
                        gpx = N*mXsCrd;
                        
                        %--------------------------------------------------
                        % State 2 (auxiliary)
                        %--------------------------------------------------
                        r = sqrt((gpx(:,1)-vTpCrd(1)).^2+...
                                 (gpx(:,2)-vTpCrd(2)).^2);
                        
                        [~,dudx] = CrackTipField_Jump(K1,K2,G_tip,k_tip,r);
                        %--------------------------------------------------
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2);
                        
                        for i = 1:ngp
                            Ms_tip = Ms_tip + dudx(i,:)*p*dsdX*w(i);
                        end
                        
                    end
                end
                
                %----------------------------------------------------------
                % Loop over boundary elements
                %----------------------------------------------------------
                
                for i_bnd = 1:length(jElBnd)
                    
                    mElCrd = mNdCrd(mLNodS(vElEnr(jElBnd(i_bnd)),:),:);
                    
                    % get intersection of crack with el.
                    [nXsCrd,mXsCrd] = GeoXElem(mSgCrd,mElCrd);
                    
                    % to integrate, must have a length and must have a jump 
                    if nXsCrd == 2 && ~all(LevelSet2Line(mElCrd,mSgCrd)>-tol_abs)
                        
                        gpx = N*mXsCrd;
                        
                        if with_MapTyp == 1
                            % local mapping
                            gpX = N*Gauss_Glb2Lcl(mXsCrd,mElCrd);
                        else
                            % global mapping
                            gpX = Gauss_Glb2Lcl(gpx,mElCrd);
                        end                        

                        q = ShapesStd_omg(ngp,gpX)*mBnWgt(:,i_bnd);
                        
                        %--------------------------------------------------
                        % State 2 (auxiliary)
                        %--------------------------------------------------
                        r = sqrt((gpx(:,1)-vTpCrd(1)).^2+...
                                 (gpx(:,2)-vTpCrd(2)).^2);
                        
                        [~,dudx] = CrackTipField_Jump(K1,K2,G_tip,k_tip,r);
                        %--------------------------------------------------
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2); % = ds/dX (2D)
                        
                        for i = 1:ngp
                            Ms_tip = Ms_tip + dudx(i,:)*p*q(i)*dsdX*w(i);
                        end
                        
                    end
                end
                
                Ms(i_crk,i_tip) = -Ms_tip;
                
            end
        end
    end
end
    
%==========================================================================
    
end