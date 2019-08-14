
%==========================================================================
% Enriched element quadrature
%==========================================================================

switch mesh_ElemType
    
    case 'T3'
        
        %------------------------------------------------------------------
        % Brn.
        %------------------------------------------------------------------
        
%         % col.
%         [mGsPnt,vGsBrn_colWgt] = Gauss_DomQuad(nGsBrn_col);
%         [mGsBrn_colShp,mGsBrn_colDrv] = ShapesStd_quad(nGsBrn_col,mGsPnt);
        
        % omg.
        [mGsPnt,vGsBrn_omgWgt] = Gauss_DomTri(nGsBrn_omg);
        [mGsBrn_omgShp,mGsBrn_omgDrv] = ShapesStd_tri(nGsBrn_omg,mGsPnt);
        
        % sub.
        [mGsPnt,vGsBrn_subWgt] = Gauss_DomTri(nGsBrn_sub);
        [mGsBrn_subShp,mGsBrn_subDrv] = ShapesStd_tri(nGsBrn_sub,mGsPnt);
        
        % gam.
        [mGsPnt,vGsBrn_gamWgt] = Gauss_DomLin(nGsBrn_gam);
        [mGsBrn_gamShp,mGsBrn_gamDrv] = ShapesStd_lin(nGsBrn_gam,mGsPnt);
        
        %------------------------------------------------------------------
        % Hvi.
        %------------------------------------------------------------------
        
        % sub
        [mGsPnt,vGsHvi_subWgt] = Gauss_DomTri(nGsHvi_sub);
        [mGsHvi_subShp,mGsHvi_subDrv] = ShapesStd_tri(nGsHvi_sub,mGsPnt);
        
        % gam
        [mGsPnt,vGsHvi_gamWgt] = Gauss_DomLin(nGsHvi_gam);
        [mGsHvi_gamShp,mGsHvi_gamDrv] = ShapesStd_lin(nGsHvi_gam,mGsPnt);
        
    case 'Q4'
        
        %------------------------------------------------------------------
        % Brn.
        %------------------------------------------------------------------
        
%         % col.
%         [mGsPnt,vGsBrn_colWgt] = Gauss_DomQuad(nGsBrn_col);
%         [mGsBrn_colShp,mGsBrn_colDrv] = ShapesStd_quad(nGsBrn_col,mGsPnt);
        
        % omg.
        [mGsPnt,vGsBrn_omgWgt] = Gauss_DomQuad(nGsBrn_omg);
        [mGsBrn_omgShp,mGsBrn_omgDrv] = ShapesStd_quad(nGsBrn_omg,mGsPnt);
        
        % sub.
        [mGsPnt,vGsBrn_subWgt] = Gauss_DomTri(nGsBrn_sub);
        [mGsBrn_subShp,mGsBrn_subDrv] = ShapesStd_tri(nGsBrn_sub,mGsPnt);
        
        % gam.
        [mGsPnt,vGsBrn_gamWgt] = Gauss_DomLin(nGsBrn_gam);
        [mGsBrn_gamShp,mGsBrn_gamDrv] = ShapesStd_lin(nGsBrn_gam,mGsPnt);
        
        %------------------------------------------------------------------
        % Hvi.
        %------------------------------------------------------------------
        
        % sub
        [mGsPnt,vGsHvi_subWgt] = Gauss_DomTri(nGsHvi_sub);
        [mGsHvi_subShp,mGsHvi_subDrv] = ShapesStd_tri(nGsHvi_sub,mGsPnt);
        
        % gam
        [mGsPnt,vGsHvi_gamWgt] = Gauss_DomLin(nGsHvi_gam);
        [mGsHvi_gamShp,mGsHvi_gamDrv] = ShapesStd_lin(nGsHvi_gam,mGsPnt);
        
    otherwise
        
        error('Unknown element type: %s',mesh_ElemType)
        
end
