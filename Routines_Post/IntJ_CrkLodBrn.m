function Gs = IntJ_CrkLodBrn(cCkCrd,mNDspE,ngp,N,w,P,R)

%==========================================================================

% Crack Tip Energy Release Rate: The J-Integral

% direction of derivative needs to point towards crack tip

%==========================================================================

%--------------------------------------------------------------------------
% Define Global Variables
%--------------------------------------------------------------------------

global mNdCrd
global mLNodS

global cLNodE
global cLEnDt

global cBrStd_eTp
global cBrStd_eSp
global cBrBln_eSp
global cBrBln_rSp

global with_MapTyp
global tol_abs

%--------------------------------------------------------------------------

jEnFun_jmp = 1;
nLNodS_enr = size(mLNodS,2);

%--------------------------------------------------------------------------
jLNodS_enr = 1:nLNodS_enr;
nEnFun_jmp = length(jEnFun_jmp);
jLNodE_zro = zeros(nEnFun_jmp,nLNodS_enr);

for iEnFun = 1:nEnFun_jmp
    jLNodE_zro(iEnFun,:) = jLNodS_enr + nLNodS_enr*(jEnFun_jmp(iEnFun)-1);
end
%--------------------------------------------------------------------------

nDimes = 2;

Gs = zeros(length(cCkCrd),2);

dXdx = zeros(2,2);
T    = zeros(2,2);

i0_gdv = 1:nDimes;

for i_crk = 1:length(cCkCrd)
    
    p = P(i_crk,[2,1])'; % [pressure,shear] -> [shear,pressure] (convenient)
    
    if p(1)^2+p(2)^2 > 0
        
        nCkCrd = size(cCkCrd{i_crk},1);
        jSgCrd = [2,1;nCkCrd-1,nCkCrd]; % (towards tip)
        
        for i_tip = 1:2
            if R(i_crk,i_tip) > 0 && ~isempty(cBrStd_eTp{i_crk,i_tip})
                
                Gs_tip = 0;
                
                mSgCrd = cCkCrd{i_crk}(jSgCrd(i_tip,:),:);
                vTpCrd = mSgCrd(2,:);
                
                % dir. of segment
                n = vTpCrd-mSgCrd(1,:);
                n = n/sqrt(n(1)^2+n(2)^2);
                
                % local to global transf.
                T = [n(1),-n(2);n(2),n(1)];
                
                t = T*p; % (global tractions)
                
                %----------------------------------------------------------
                % Loop over std. elements
                %----------------------------------------------------------
                
                % get candidate elements
                vElEnr = [cBrStd_eTp{i_crk,i_tip};cBrStd_eSp{i_crk,i_tip}];
                
                % detect elements within domain; get smoothing weight
                [jElDom,jElBnd,mBnWgt] = ...
                    Elems_tip(mNdCrd,mLNodS(vElEnr,:),vTpCrd,R(i_crk,i_tip));
                
                % good to transpose here
                mBnWgt = mBnWgt';
                
                for i_dom = jElDom(:)'
                    i_elm = vElEnr(i_dom);
                    
                    % el. enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % enr. type for Brn. == i_tip
                    i = find(mLEnDt(1,:) == i_crk);
                    i = i(mLEnDt(2,i) == i_tip);
                    
                    % get 0th enr. node
                    iLNodE = sum(mLEnDt(3,1:i(1)-1)); % sum([]) == 0 (okey)
                    
                    mElCrd = mNdCrd(mLNodS(i_elm,:),:);
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
                        
                        [mGsShp,mGsDrv] = ShapesStd_omg(ngp,gpX);
                        [f,dfds] = EnrJmp_Brn(gpx,vTpCrd);
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2);
                        dxdX = mGsDrv*mElCrd;

                        for i_enr = 1:nEnFun_jmp
                            
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mLDspE = mNDspE(:,cLNodE{i_elm}(jLNodE))';
                            
                            u    = mGsShp*mLDspE;
                            dudX = mGsDrv*mLDspE; % later multiply by jump.
                            
                            i_gdv = i0_gdv;
                            for i = 1:ngp
                                
                                J       = dxdX(i_gdv,:);
                                detJ    = J(1)*J(4)-J(2)*J(3);
                                
                                dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                                dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                                
                                ux = f(i,i_enr)*n*dXdx*dudX(i_gdv,:) + ...
                                     dfds(i,i_enr)*u(i,:);
                                 
                                Gs_tip = Gs_tip + ux*t*dsdX*w(i);
                                
                                i_gdv = i_gdv + nDimes;
                                
                            end
                        end
                    end
                end
                
                for i_bnd = 1:length(jElBnd)
                    i_elm = vElEnr(jElBnd(i_bnd));
                    
                    % el. enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % enr. type for Brn. == i_tip
                    i = find(mLEnDt(1,:) == i_crk);
                    i = i(mLEnDt(2,i) == i_tip);
                    
                    % get 0th enr. node
                    iLNodE = sum(mLEnDt(3,1:i(1)-1)); % sum([]) == 0 (okey)
                    
                    mElCrd = mNdCrd(mLNodS(i_elm,:),:);
                    [nXsCrd,mXsCrd] = GeoXElem(mSgCrd,mElCrd);
                    
                    if nXsCrd == 2 && ~all(LevelSet2Line(mElCrd,mSgCrd)>-tol_abs)
                        
                        gpx = N*mXsCrd;
                        
                        if with_MapTyp == 1
                            % local mapping
                            gpX = N*Gauss_Glb2Lcl(mXsCrd,mElCrd);
                        else
                            % global mapping
                            gpX = Gauss_Glb2Lcl(gpx,mElCrd);
                        end
                        
                        [mGsShp,mGsDrv] = ShapesStd_omg(ngp,gpX);
                        [f,dfds] = EnrJmp_Brn(gpx,vTpCrd);
                        
                        q = mGsShp*mBnWgt(:,i_bnd);
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2);
                        dxdX = mGsDrv*mElCrd;
                        
                        for i_enr = 1:nEnFun_jmp
                        
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mLDspE = mNDspE(:,cLNodE{i_elm}(jLNodE))';
                            
                            u    = mGsShp*mLDspE;
                            dudX = mGsDrv*mLDspE; % later multiply by jump.
                            
                            i_gdv = i0_gdv;
                            for i = 1:ngp
                                
                                J       = dxdX(i_gdv,:);
                                detJ    = J(1)*J(4)-J(2)*J(3);
                                
                                dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                                dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                                
                                ux = f(i,i_enr)*n*dXdx*dudX(i_gdv,:) + ...
                                     dfds(i,i_enr)*u(i,:);
                                
                                Gs_tip = Gs_tip + ux*t*q(i)*dsdX*w(i);
                                
                                i_gdv = i_gdv + nDimes;
                                
                            end
                        end
                    end
                end
                
                %----------------------------------------------------------
                % Loop over bln. elements
                %----------------------------------------------------------
                
                % get candidate elements and blending weights
                vElEnr = cBrBln_eSp{i_crk,i_tip};
                mBlWgt = cBrBln_rSp{i_crk,i_tip}';
                
                % detect elements within domain; get smoothing weight
                [jElDom,jElBnd,mBnWgt] = ...
                    Elems_tip(mNdCrd,mLNodS(vElEnr,:),vTpCrd,R(i_crk,i_tip));
                
                % good to transpose here
                mBnWgt = mBnWgt';
                
                for i_dom = jElDom(:)'
                    i_elm = vElEnr(i_dom);
                    
                    % el. enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % enr. type for Brn. == i_tip
                    i = find(mLEnDt(1,:) == i_crk);
                    i = i(mLEnDt(2,i) == i_tip);
                    
                    % get 0th enr. node
                    iLNodE = sum(mLEnDt(3,1:i(1)-1)); % sum([]) == 0 (okey)
                    
                    mElCrd = mNdCrd(mLNodS(i_elm,:),:);
                    [nXsCrd,mXsCrd] = GeoXElem(mSgCrd,mElCrd);
                    
                    if nXsCrd == 2 && ~all(LevelSet2Line(mElCrd,mSgCrd)>-tol_abs)
                        
                        gpx = N*mXsCrd;
                        
                        if with_MapTyp == 1
                            % local mapping
                            gpX = N*Gauss_Glb2Lcl(mXsCrd,mElCrd);
                        else
                            % global mapping
                            gpX = Gauss_Glb2Lcl(gpx,mElCrd);
                        end
                        
                        [mGsShp,mGsDrv] = ShapesStd_omg(ngp,gpX);
                        [f,dfds] = EnrJmp_Brn(gpx,vTpCrd);
                        
                        W    = mGsShp*mBlWgt(:,i_dom);
                        dWdX = mGsDrv*mBlWgt(:,i_dom);
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2);
                        dxdX = mGsDrv*mElCrd;
                        
                        for i_enr = 1:nEnFun_jmp
                            
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mLDspE = mNDspE(:,cLNodE{i_elm}(jLNodE))';
                            
                            u    = mGsShp*mLDspE;
                            dudX = mGsDrv*mLDspE; % later multiply by jump.
                            
                            i_gdv = i0_gdv;
                            for i = 1:ngp
                                
                                J       = dxdX(i_gdv,:);
                                detJ    = J(1)*J(4)-J(2)*J(3);
                                
                                dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                                dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                                
                                ux = f(i,i_enr)*W(i)*n*dXdx*dudX(i_gdv,:) + ...
                                     f(i,i_enr)*n*dXdx*dWdX(i_gdv)*u(i,:) + ...
                                     dfds(i,i_enr)*W(i)*u(i,:);
                                 
                                Gs_tip = Gs_tip + ux*t*dsdX*w(i);
                                
                                i_gdv = i_gdv + nDimes;
                                
                            end
                            
                        end
                    end
                end
                
                for i2_bnd = 1:length(jElBnd)
                    i1_bnd = jElBnd(i2_bnd);
                    i_elm = vElEnr(i1_bnd);
                    
                    % el. enr. data
                    mLEnDt = cLEnDt{i_elm};
                    
                    % enr. type for Brn. == i_tip
                    i = find(mLEnDt(1,:) == i_crk);
                    i = i(mLEnDt(2,i) == i_tip);
                    
                    % get 0th enr. node
                    iLNodE = sum(mLEnDt(3,1:i(1)-1)); % sum([]) == 0 (okey)
                    
                    mElCrd = mNdCrd(mLNodS(i_elm,:),:);
                    [nXsCrd,mXsCrd] = GeoXElem(mSgCrd,mElCrd);
                    
                    if nXsCrd == 2 && ~all(LevelSet2Line(mElCrd,mSgCrd)>-tol_abs)
                        
                        gpx = N*mXsCrd;
                        
                        if with_MapTyp == 1
                            % local mapping
                            gpX = N*Gauss_Glb2Lcl(mXsCrd,mElCrd);
                        else
                            % global mapping
                            gpX = Gauss_Glb2Lcl(gpx,mElCrd);
                        end
                        
                        [mGsShp,mGsDrv] = ShapesStd_omg(ngp,gpX);
                        [f,dfds] = EnrJmp_Brn(gpx,vTpCrd);
                        
                        q    = mGsShp*mBnWgt(:,i2_bnd);
                        W    = mGsShp*mBlWgt(:,i1_bnd);
                        dWdX = mGsDrv*mBlWgt(:,i1_bnd);
                        
                        d = mXsCrd(2,:)-mXsCrd(1,:);
                        dsdX = 0.5*sqrt(d(1)^2+d(2)^2);
                        dxdX = mGsDrv*mElCrd;
                        
                        for i_enr = 1:nEnFun_jmp
                            
                            jLNodE = iLNodE + jLNodE_zro(i_enr,:);
                            mLDspE = mNDspE(:,cLNodE{i_elm}(jLNodE))';
                            
                            u    = mGsShp*mLDspE;
                            dudX = mGsDrv*mLDspE; % later multiply by jump.
                            
                            i_gdv = i0_gdv;
                            for i = 1:ngp
                                
                                J       = dxdX(i_gdv,:);
                                detJ    = J(1)*J(4)-J(2)*J(3);
                                
                                dXdx(1) =  J(4)/detJ; dXdx(3) = -J(3)/detJ;
                                dXdx(2) = -J(2)/detJ; dXdx(4) =  J(1)/detJ;
                                
                                ux = f(i,i_enr)*W(i)*n*dXdx*dudX(i_gdv,:) + ...
                                     f(i,i_enr)*n*dXdx*dWdX(i_gdv)*u(i,:) + ...
                                     dfds(i,i_enr)*W(i)*u(i,:);
                                 
                                Gs_tip = Gs_tip + ux*t*q(i)*dsdX*w(i);
                                
                                i_gdv = i_gdv + nDimes;
                                
                            end
                        end
                    end
                end
                
                Gs(i_crk,i_tip) = -Gs_tip;
                
            end
        end
    end
end
    
%==========================================================================
    
end