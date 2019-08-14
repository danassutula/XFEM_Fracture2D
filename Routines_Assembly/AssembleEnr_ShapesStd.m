
%==========================================================================
% Determine standard shape functions for enriched elements
%
%   N.B.:   need to follow in order of presedence
%           e.g. crack tip quadrature overrides other types
%           e.g. Branch quadrature overrides Heaviside quadrature
%
%==========================================================================


%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

cGsEnr_omgShS = cell(nElems,1);
cGsEnr_omgDvS = cell(nElems,1);
cGsEnr_omgWgt = cell(nElems,1);

%--------------------------------------------------------------------------
% Branch enrichment: crack tip elements
%--------------------------------------------------------------------------

for i = 1:nCrack
for j = 1:2
for k = cBrStd_eTp{i,j}';
if  isempty(cGsEnr_omgWgt{k}); mElCrd=mNdCrd(mLNodS(k,:),:);
    
    if length(Nodes_unq(cXsElm{k,1},tol_abs)) < 3
        
        if with_MapTyp == 1
        
        %------------------------------------------------------------------
        % Gauss points and weights (sub-div. local)
        %------------------------------------------------------------------
        
        mSbCrd = [Gauss_Glb2Lcl(cXsElm{k,1},mElCrd);mElCrd_loc];
        mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_ref),:);
        mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip
        
        % (op1) full-polar
        [nGauss,mGsPnt,vGsWgt] = ...
            Gauss_SubDivLcl_Plr(mSbCrd,mSbNod,nGsBrn_pol);
        
%         % (op2) quasi-polar
%         mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%         [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
%            nGsBrn_col,mGsBrn_colShp,mGsBrn_colDrv,vGsBrn_colWgt);
        
        else
        
        %------------------------------------------------------------------
        % Gauss points and weights (sub-div. global)
        %------------------------------------------------------------------
        
        mSbCrd = [cXsElm{k,1};mElCrd];
        mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_abs),:);
        mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip
        
        % (op1) full-polar
        [nGauss,mGsPnt,vGsWgt] = ...
            Gauss_SubDivGlb_Plr(mElCrd,mSbCrd,mSbNod,nGsBrn_pol);
        
%         % (op2) quasi-polar
%         mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%         [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
%             nGsBrn_col,mGsBrn_colShp,mGsBrn_colDrv,vGsBrn_colWgt);
        
        end
        
    else % tip element subdivided by multiple cracks
        
        if with_MapTyp == 1
            
        %------------------------------------------------------------------
        % Gauss points and weights (sub-div. local)
        %------------------------------------------------------------------
        
        mSbCrd = [Gauss_Glb2Lcl(cXsElm{k,1},mElCrd);mElCrd_loc];
        [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{k,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
            nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
        
        else
            
        %------------------------------------------------------------------
        % Gauss points and weights (sub-div. global)
        %------------------------------------------------------------------
        
        [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{k,1};mElCrd],cXsElm{k,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
            nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
        
        end
        
    end
    
    %----------------------------------------------------------------------
    % Gauss shapes (std.)
    %----------------------------------------------------------------------
    
    [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
    
    cGsEnr_omgShS{k} = mGsShp;
    cGsEnr_omgDvS{k} = mGsDrv;
    cGsEnr_omgWgt{k} = vGsWgt;
    
end
end
end
end

%--------------------------------------------------------------------------
% Branch enrichment: split standard & blending elements (can be grouped)
%--------------------------------------------------------------------------

for i = 1:nCrack
for j = 1:2
for k = [cBrStd_eSp{i,j};cBrBln_eSp{i,j}]'
if  isempty(cGsEnr_omgWgt{k}); mElCrd=mNdCrd(mLNodS(k,:),:);
    
    if with_MapTyp == 1
    
    mSbCrd = [Gauss_Glb2Lcl(cXsElm{k,1},mElCrd);mElCrd_loc];
    [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{k,2});
    [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
        nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
    
    else
    
    [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{k,1};mElCrd],cXsElm{k,2});
    [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
        nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
    
    end
    
    [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
    
    cGsEnr_omgShS{k} = mGsShp;
    cGsEnr_omgDvS{k} = mGsDrv;
    cGsEnr_omgWgt{k} = vGsWgt;
    
end
end
end
end


%--------------------------------------------------------------------------
% Branch enrichment: full standard & blending elements (can be grouped)
%--------------------------------------------------------------------------

for i = 1:nCrack
for j = 1:2
    
    if j == 1
        x_tip = cCkCrd{i}(1,:);
    else % j == 2
        x_tip = cCkCrd{i}(end,:);
    end
    
for k = [cBrStd_eFl{i,j};cBrBln_eFl{i,j}]'
if  isempty(cGsEnr_omgWgt{k}); mElCrd=mNdCrd(mLNodS(k,:),:);
    
    if isempty(cXsElm{k,1})
        
        [x_lvl,s_lvl] = LevelSet2Elem(x_tip,mElCrd);
        if s_lvl < 0.1*ElemSize(mElCrd) % tip too close to next el.

            if with_MapTyp == 1

            mSbCrd = [Gauss_Glb2Lcl(x_lvl,mElCrd); mElCrd_loc];
            mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_ref),:);
            mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip

            % (op1) full-polar
            [nGauss,mGsPnt,vGsWgt] = ...
                Gauss_SubDivLcl_Plr(mSbCrd,mSbNod,nGsBrn_pol);

%             % (op2) quasi-polar
%             mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%             [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
%                 nGsBrn_omg,mGsBrn_omgShp,mGsBrn_omgDrv,vGsBrn_omgWgt);

            else

            mSbCrd = [x_lvl;mElCrd];
            mSbCrd = mSbCrd(Nodes_unq(mSbCrd,tol_abs),:);
            mSbNod = GeoSubTri_Center(mSbCrd); % 1st node is crack tip
            
            % (op1) full-polar
            [nGauss,mGsPnt,vGsWgt] = ...
                Gauss_SubDivGlb_Plr(mElCrd,mSbCrd,mSbNod,nGsBrn_pol);
            
%             % (op2) quasi-polar
%             mSbNod(:,4) = mSbNod(:,1); % convert tri. to quad.
%             [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mSbCrd,mSbNod,...
%                 nGsBrn_omg,mGsBrn_omgShp,mGsBrn_omgDrv,vGsBrn_omgWgt);
            
            end
            
            [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
            
            cGsEnr_omgShS{k} = mGsShp;
            cGsEnr_omgDvS{k} = mGsDrv;
            cGsEnr_omgWgt{k} = vGsWgt;
            
        else
            
            cGsEnr_omgShS{k} = mGsBrn_omgShp;
            cGsEnr_omgDvS{k} = mGsBrn_omgDrv;
            cGsEnr_omgWgt{k} = vGsBrn_omgWgt;
            
        end
        
    else

        if with_MapTyp == 1
        
        mSbCrd = [Gauss_Glb2Lcl(cXsElm{k,1},mElCrd);mElCrd_loc];
        [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{k,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
            nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
        
        else
        
        [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{k,1};mElCrd],cXsElm{k,2});
        [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
            nGsBrn_sub,mGsBrn_subShp,mGsBrn_subDrv,vGsBrn_subWgt);
        
        end
        
        [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
        
        cGsEnr_omgShS{k} = mGsShp;
        cGsEnr_omgDvS{k} = mGsDrv;
        cGsEnr_omgWgt{k} = vGsWgt;
        
    end
    
end
end
end
end

%--------------------------------------------------------------------------
% Heaviside enrichment: standard & blending elements (can be grouped)
%--------------------------------------------------------------------------

for i = 1:nCrack
for k = [cHvBln_elm{i,1};cHvStd_elm{i};cHvBln_elm{i,2}]'
if  isempty(cGsEnr_omgWgt{k}); mElCrd=mNdCrd(mLNodS(k,:),:);
    
    if with_MapTyp == 1
    
    mSbCrd = [Gauss_Glb2Lcl(cXsElm{k,1},mElCrd);mElCrd_loc];
    [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,cXsElm{k,2});
    [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivLcl(mSbCrd,mSbNod,...
        nGsHvi_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt);
    
    else
    
    [mSbCrd,mSbNod] = GeoSubTri_Delaun([cXsElm{k,1};mElCrd],cXsElm{k,2});
    [nGauss,mGsPnt,vGsWgt] = Gauss_SubDivGlb(mElCrd,mSbCrd,mSbNod,...
        nGsHvi_sub,mGsHvi_subShp,mGsHvi_subDrv,vGsHvi_subWgt);
    
    end
    
    [mGsShp,mGsDrv] = ShapesStd_omg(nGauss,mGsPnt);
    
    cGsEnr_omgShS{k} = mGsShp;
    cGsEnr_omgDvS{k} = mGsDrv;
    cGsEnr_omgWgt{k} = vGsWgt;
    
end
end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% OVERWRITE ALL QUADRATURE USING A HIGH ORDER GAUSS RULE
% - an equivalent script appears in 'AssembleUpd_ShapesStd'
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% NOT ROBUST. THE STIFFNESS MATRIX IS SOMETIMES SINGULAR.

if 0
    
    % use a very high order quadrature
    
    q = false(nElEnr,1);
    for k = vElEnr(:)'
        if length(Nodes_unq(cXsElm{k},tol_abs))>4
            % at least two cracks cut the element
            q(vElEnr==k) = 1;
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
        
        cGsEnr_omgShS(vElEnr(q)) = {mGsShp};
        cGsEnr_omgDvS(vElEnr(q)) = {mGsDrv};
        cGsEnr_omgWgt(vElEnr(q)) = {vGsWgt};
        
    end
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
