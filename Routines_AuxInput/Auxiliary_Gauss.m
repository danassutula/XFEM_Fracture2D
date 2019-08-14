
%--------------------------------------------------------------------------
% Gauss Quadrature
%--------------------------------------------------------------------------

switch mesh_ElemType

    case 'T3'
        
        % std.
        nGsStd_omg = 1;
        nGsStd_gam = 1;
        
        % brn.: tip element(s)
        nGsBrn_pol = [13,7];
        
        % brn.: non-split elements (std. & bln.)
        nGsBrn_omg = 33; % 19; % 33
        
        % brn.: split elements (std. & bln.)
        nGsBrn_sub = 33; % 19 % 33
        
        % brn.: surface boundary
        nGsBrn_gam = 49;
        
        % hvi.: split elements (std. & bln.)
        nGsHvi_sub = 1;
        
        % hvi.: surface boundary
        nGsHvi_gam = 1;
        
    case 'Q4'
        
        % std
        nGsStd_omg = 4;
        nGsStd_gam = 2;
        
        % brn.: tip element(s)
        nGsBrn_pol = [13,7];
        
        % brn.: non-split elements (std. & bln.)
        nGsBrn_omg = 36; % 49; 36; 25
        
        % brn.: split elements (std. & bln.)
        nGsBrn_sub = 33; % 19
        
        % brn.: surface boundary
        nGsBrn_gam = 49;
        
        % hvi.: split elements (std. & bln.)
        nGsHvi_sub = 3; % 33; % 3;
        
        % hvi.: surface boundary
        nGsHvi_gam = 2;
        
    otherwise
        
        error('Unknown element type: %s',mesh_ElemType)

end

%--------------------------------------------------------------------------
