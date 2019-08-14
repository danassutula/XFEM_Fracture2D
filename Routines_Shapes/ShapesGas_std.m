
%==========================================================================
% Standard element quadrature
%==========================================================================

switch mesh_ElemType
    
    case 'T3'
        
        % omg.
        [mGsPnt,vGsStd_omgWgt] = Gauss_DomTri(nGsStd_omg);
        [mGsStd_omgShp,mGsStd_omgDrv] = ShapesStd_tri(nGsStd_omg,mGsPnt);
        
        % gam.
        [mGsPnt,vGsStd_gamWgt] = Gauss_DomLin(nGsStd_gam);
        [mGsStd_gamShp,mGsStd_gamDrv] = ShapesStd_lin(nGsStd_gam,mGsPnt);
        
    case 'Q4'
        
        % omg.
        [mGsPnt,vGsStd_omgWgt] = Gauss_DomQuad(nGsStd_omg);
        [mGsStd_omgShp,mGsStd_omgDrv] = ShapesStd_quad(nGsStd_omg,mGsPnt);
        
        % gam.
        [mGsPnt,vGsStd_gamWgt] = Gauss_DomLin(nGsStd_gam);
        [mGsStd_gamShp,mGsStd_gamDrv] = ShapesStd_lin(nGsStd_gam,mGsPnt);
        
    otherwise
        
        error('Unknown element type: %s',mesh_ElemenType)
        
end