
%==========================================================================
% Determine standard shape functions for enriched elements
%==========================================================================

%--------------------------------------------------------------------------
% Initialize (already initialized)
%--------------------------------------------------------------------------
% cGsEnr_omgShS = cell(nElems,1);
% cGsEnr_omgDvS = cell(nElems,1);
% cGsEnr_omgWgt = cell(nElems,1);
%--------------------------------------------------------------------------

for i = vElUpd(:)'
    
    mElCrd = mNdCrd(mLNodS(i,:),:);
    vLQTyp = cLEnDt{i}(4,:);
    
switch max(vLQTyp)
    case 1 % a Heaviside element: std. & bln.
        
        if with_MapTyp == 1
        
        mSbCrd = [Gauss_Glb2Lcl(cXsElm{i,1},mElCrd);mElCrd_loc];
        [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{i,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
            nGsHvi_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt);
        
        else
        
        [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{i,1};mElCrd],cXsElm{i,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
            nGsHvi_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt);
        
        end
        
        [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
        
        cGsEnr_omgShS{i} = mGsShp;
        cGsEnr_omgDvS{i} = mGsDrv;
        cGsEnr_omgWgt{i} = vGsWgt;
        
    case 2 % a Branch element: std. & bln. (no crack tip elements)
        
        if isempty(cXsElm{i})
            
            if cLEnDt{i}(2) == 1
                x_tip = cCkCrd{cLEnDt{i}(1)}(1,:);
            else % cLEnDt{i}(2) == 2
                x_tip = cCkCrd{cLEnDt{i}(1)}(end,:);
            end
            
            [x_lvl,s_lvl] = LevelSet2Elem(x_tip,mElCrd);
            if s_lvl < 0.1*ElemSize(mElCrd) % tip too close to next el.
            
                if with_MapTyp == 1
                
                mSbCrd = [Gauss_Glb2Lcl(x_lvl,mElCrd);mElCrd_loc];
                mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_ref),:);
                mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip

                % (op1) full-polar
                [nGauss,mGsPnt,vGsWgt] = ...
                    Gauss_SubDivLcl_Plr(mSbCrd,mSbNod,nGsBrn_pol);

%                 % (op2) quasi-polar
%                 mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%                 [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
%                     nGsBrn_omg,mGsBrn_omgShp,mGsBrn_omgDrv,vGsBrn_omgWgt);

                else
                
                mSbCrd = [x_lvl;mElCrd];
                mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_abs),:);
                mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip

                % (op1) full-polar
                [nGauss,mGsPnt,vGsWgt] = ...
                    Gauss_SubDivGlb_Plr(mElCrd,mSbCrd,mSbNod,nGsBrn_pol);
                
%                 % (op2) quasi-polar
%                 mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%                 [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mSbCrd,mSbNod,...
%                     nGsBrn_omg,mGsBrn_omgShp,mGsBrn_omgDrv,vGsBrn_omgWgt);
                    
                end
                
                [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);

                cGsEnr_omgShS{i} = mGsShp;
                cGsEnr_omgDvS{i} = mGsDrv;
                cGsEnr_omgWgt{i} = vGsWgt;
                
            else
                
                cGsEnr_omgShS{i} = mGsBrn_omgShp;
                cGsEnr_omgDvS{i} = mGsBrn_omgDrv;
                cGsEnr_omgWgt{i} = vGsBrn_omgWgt;
            
            end
            
        else
            
            if with_MapTyp == 1
            
            mSbCrd = [Gauss_Glb2Lcl(cXsElm{i,1},mElCrd);mElCrd_loc];
            [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{i,2});
            [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
                nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
            
            else
            
            [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{i,1};mElCrd],cXsElm{i,2});
            [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
                nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
            
            end
            
            [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
            
            cGsEnr_omgShS{i} = mGsShp;
            cGsEnr_omgDvS{i} = mGsDrv;
            cGsEnr_omgWgt{i} = vGsWgt;
            
        end
        
    case 3 % a Branch element: crack tip element
        
        if length(Nodes_unq(cXsElm{i},tol_abs)) < 3
            
            if with_MapTyp == 1
            
            mSbCrd = [Gauss_Glb2Lcl(cXsElm{i,1},mElCrd);mElCrd_loc];
            mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_ref),:);
            mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip

            % (op1) full-polar
            [nGauss,mGsPnt,vGsWgt] = ...
                Gauss_SubDivLcl_Plr(mSbCrd,mSbNod,nGsBrn_pol);

%             % (op2) quasi-polar
%             mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%             [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
%                nGsBrn_col,mGsBrn_colShp,mGsBrn_colDrv,vGsBrn_colWgt);
            
            else
            
            mSbCrd = [cXsElm{i,1};mElCrd];
            mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_abs),:);
            mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip

            % (op1) full-polar
            [nGauss,mGsPnt,vGsWgt] = ...
                Gauss_SubDivGlb_Plr(mElCrd,mSbCrd,mSbNod,nGsBrn_pol);
            
%             % (op2) quasi-polar
%             mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%             [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
%                 nGsBrn_col,mGsBrn_colShp,mGsBrn_colDrv,vGsBrn_colWgt);
            
            end
            
        else % integration over sub-cells
            
            if with_MapTyp == 1
            
            mSbCrd = [Gauss_Glb2Lcl(cXsElm{i,1},mElCrd);mElCrd_loc];
            [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{i,2});
            [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
                nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
            
            else
            
            [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{i,1};mElCrd],cXsElm{i,2});
            [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
                nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
            
            end
            
        end
        
        [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
        
        cGsEnr_omgShS{i} = mGsShp;
        cGsEnr_omgDvS{i} = mGsDrv;
        cGsEnr_omgWgt{i} = vGsWgt;
        
    otherwise
        error('unknown integration rule for element: %d',i)
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% OVERWRITE ALL QUADRATURE USING A HIGH ORDER GAUSS RULE
% - an equivalent script appears in 'AssembleEnr_ShapesStd'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% NOT ROBUST. THE STIFFNESS MATRIX IS SOMETIMES SINGULAR.

if 0
    
    % use a very high order quadrature
    
    q = false(nElUpd,1);
    for k = vElUpd(:)'
        if length(Nodes_unq(cXsElm{k},tol_abs))>4
            % at least two cracks cut the element
            q(vElUpd==k) = 1;
        end
    end
    
    if any(q)
        
        warning('using a very high order quadrature for multiply-intersected element(s)')
        
        switch mesh_ElemType
            case 'T3'; nGauss = 33;
                [mGsPnt,vGsWgt] = Gauss_DomTri(nGauss);
                [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
                
            case 'Q4'; nGauss = 31^2;
                [mGsPnt,vGsWgt] = Gauss_DomQuad(nGauss);
                [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
        end
        
        cGsEnr_omgShS(vElUpd(q)) = {mGsShp};
        cGsEnr_omgDvS(vElUpd(q)) = {mGsDrv};
        cGsEnr_omgWgt(vElUpd(q)) = {vGsWgt};
        
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

