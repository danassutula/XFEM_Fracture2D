function vGFrcE = ForceEnr_BndCkBrn(nGlDof,nElEnr,jEnFun_jmp,nLNodS_enr,...
    nCrack,cCkCrd,nGauss,mGsShp_gam,vGsWgt,mCkLod)

% vCkLod(iCrack,:) = [ p , \tau ]

%--------------------------------------------------------------------------
% Define Global Variables
%--------------------------------------------------------------------------

global mNdCrd
global mLNodS

global cLEnDt
global cLNodE
global mNDofE

global cBrStd_eTp
global cBrStd_eSp
global cBrBln_eSp
global cBrBln_rSp

global with_MapTyp
global tol_abs

%--------------------------------------------------------------------------

jLNodS_enr = 1:nLNodS_enr;
nEnFun_jmp = length(jEnFun_jmp);
jLNodE_zro = zeros(nEnFun_jmp,nLNodS_enr);

for i_enr = 1:nEnFun_jmp
    jLNodE_zro(i_enr,:) = jLNodS_enr + nLNodS_enr*(jEnFun_jmp(i_enr)-1);
end

nDimes = 2;

nElDof = nDimes*nLNodS_enr;
nSprGF = nElEnr*nEnFun_jmp*nElDof; % predictor (over estimate)
jSprGF = zeros(nSprGF,1);
vSprGF = zeros(nSprGF,1);

jSprLF = 1-nElDof:0;
jDofLF = 0:nDimes:(nElDof-nDimes);
vLFrcE = zeros(nElDof,1);

jGauss = 1:nGauss;
vShpJW_zrs = zeros(1,nLNodS_enr);

for i_crk = 1:nCrack
    p_crk = mCkLod(i_crk,:);
    
    if any(p_crk)
        
        nCkCrd = size(cCkCrd{i_crk},1);
        jSgCrd = [2,1;nCkCrd-1,nCkCrd]; % (towards tip)
        
        for i_tip = 1:2
            if ~isempty(cBrStd_eTp{i_crk,i_tip})
                
                x_sgm = cCkCrd{i_crk}(jSgCrd(i_tip,:),:);
                x_tip = x_sgm(2,:);
                
                n = x_tip-x_sgm(1,:);
                n = norm(n)\n;
                
                fx = n(1)*p_crk(2) - n(2)*p_crk(1);
                fy = n(1)*p_crk(1) + n(2)*p_crk(2);
                
                %----------------------------------------------------------
                % Standard Branch elements
                %----------------------------------------------------------
                
                vElEnr = [cBrStd_eTp{i_crk,i_tip};...
                          cBrStd_eSp{i_crk,i_tip}];
                
                for i_elE = 1:length(vElEnr)
                    i_elm = vElEnr(i_elE);
                    
                    % element enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % Enr. type for Brn. == i_tip (good coincidence)
                    tmp = find(mLEnDt(1,:) == i_crk);
                    tmp = tmp(mLEnDt(2,tmp) == i_tip); 
                    iLNodE = sum(mLEnDt(3,1:tmp(1)-1));
                              
                    x_elm = mNdCrd(mLNodS(i_elm,:),:);
                    [n_crs,x_crs] = GeoXElem(x_sgm,x_elm);
                    
                    % to integrate, must have a length and must have a jump 
                    if n_crs==2 && ~all(LevelSet2Line(x_elm,x_sgm)>-tol_abs)
                        
                        x_gsp = mGsShp_gam*x_crs;
                        
                        if with_MapTyp == 1
                            % local mapping
                            X_gsp = mGsShp_gam*Gauss_Glb2Lcl(x_crs,x_elm);
                        else
                            % global mapping
                            X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
                        end
                        
                        mGsShp = ShapesStd_omg(nGauss,X_gsp);
                        mGsJmp = EnrJmp_Brn(x_gsp,x_tip);
                        
                        s = sqrt((x_crs(2,1)-x_crs(1,1))^2 ...
                            + (x_crs(2,2)-x_crs(1,2))^2);
                        
                        detJ = 0.5*s; % = ds/dze (2D)
                        
                        for i_enr = 1:nEnFun_jmp
                            
                            vShpJW = vShpJW_zrs;
                            
                            for i_gsp = jGauss
                                vShpJW = vShpJW + ...
                                    (mGsJmp(i_gsp,i_enr)*detJ*vGsWgt(i_gsp))*mGsShp(i_gsp,:);
                            end
                            
                            vLFrcE(jDofLF+1) = vShpJW*fx;
                            vLFrcE(jDofLF+2) = vShpJW*fy;
                            
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mElDof = mNDofE(:,cLNodE{i_elm}(jLNodE));
                            
                            jSprLF = jSprLF + nElDof;
                            jSprGF(jSprLF) = mElDof(:);
                            vSprGF(jSprLF) = vLFrcE;
                            
                        end
                    end
                end
                
                %----------------------------------------------------------
                % Blending Branch elements
                %----------------------------------------------------------
                
                vElEnr = cBrBln_eSp{i_crk,i_tip};
                mBlRmp = cBrBln_rSp{i_crk,i_tip}'; % good to transpose here
                
                for i_elE = 1:length(vElEnr)
                    i_elm = vElEnr(i_elE);
                    
                    % element enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % Enr. type for Brn. == i_tip (good coincidence)
                    tmp = find(mLEnDt(1,:) == i_crk);
                    tmp = tmp(mLEnDt(2,tmp) == i_tip); 
                    iLNodE = sum(mLEnDt(3,1:tmp(1)-1));
                    
                    x_elm = mNdCrd(mLNodS(i_elm,:),:); % for geometry
                    [n_crs,x_crs] = GeoXElem(x_sgm,x_elm);
                    
                    % to integrate, must have a length and must have a jump 
                    if n_crs == 2 && ~all(LevelSet2Line(x_elm,x_sgm)>-tol_abs)
                        
                        x_gsp = mGsShp_gam*x_crs;
                        
                        if with_MapTyp == 1
                            % local mapping
                            X_gsp = mGsShp_gam*Gauss_Glb2Lcl(x_crs,x_elm);
                        else
                            % global mapping
                            X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
                        end
                        
                        mGsShp = ShapesStd_omg(nGauss,X_gsp);
                        mGsJmp = EnrJmp_Brn(x_gsp,x_tip);
                        vGsRmp = mGsShp*mBlRmp(:,i_elE);
                        
                        s = sqrt((x_crs(2,1)-x_crs(1,1))^2 ...
                            + (x_crs(2,2)-x_crs(1,2))^2);
                        
                        detJ = 0.5*s; % = ds/dze (2D)
                        
                        for i_enr = 1:nEnFun_jmp
                            
                            vShpJW = vShpJW_zrs;
                            vGsJmp = mGsJmp(:,i_enr).*vGsRmp; % weighted vGsJmp
                            
                            for i_gsp = jGauss
                                vShpJW = vShpJW + ...
                                    (vGsJmp(i_gsp)*detJ*vGsWgt(i_gsp))*mGsShp(i_gsp,:);
                            end
                            
                            vLFrcE(jDofLF+1) = vShpJW*fx;
                            vLFrcE(jDofLF+2) = vShpJW*fy;
                            
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mElDof = mNDofE(:,cLNodE{i_elm}(jLNodE));
                            
                            jSprLF = jSprLF + nElDof;
                            jSprGF(jSprLF) = mElDof(:);
                            vSprGF(jSprLF) = vLFrcE;
                            
                        end
                    end
                end
            end
        end
    end
end

nSprGF = jSprLF(end);
jSprGF = jSprGF(1:nSprGF);
vSprGF = vSprGF(1:nSprGF);

vGFrcE = sparse(jSprGF,ones(nSprGF,1),vSprGF,nGlDof,1);

end
