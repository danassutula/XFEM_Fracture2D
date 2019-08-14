
%--------------------------------------------------------------------------
% Fool-proofing input
%--------------------------------------------------------------------------

% use a different crack extension length:
% e.g. dl = max(f_inc*r_tip,with_CkIncMax)
if ~exist('with_CkIncMax','var') || isempty(with_CkIncMax)
    with_CkIncMax = 0; % do not modify crack extension lengths
end

if nCrack > length(cCkCrd)
    error('Undefined crack(s); n_crk = %d',nCrack)
end

if ~isempty(mTpFrz) && nCrack < max(mTpFrz(:,1))
    error('Undefined crack(s) to be frozen; n_crk = %d',nCrack)
end

% keep in column
cCkCrd = cCkCrd(1:nCrack);
cCkCrd = cCkCrd(:);

if nStep > 0 && nCrack > 0
    with_Growth = 1;
else
    with_Growth = 0;
    nStep = 1;
end

% assert global growth law
if with_GLwInc && ~(strcmpi(kind_GrwCrt,'all') || strcmpi(kind_GrwCrt,'custom'))
    kind_GrwCrt = 'all'; % (only a trial crack growth criterion)
end

%--------------------------------------------------------------------------
% Crack pre-intersection increment/tip refinement handling
%--------------------------------------------------------------------------

if with_RfnInc && ~with_RfnXrs
    with_RfnXrs = 1; % (override)
    nRefine_xrs = nRefine_inc;
end

% max n. bln. segments to add
if with_RfnXrs % || with_RfnInc
    nSgBln = inf; % (is safe)
else
    nSgBln = 3; % (is enough)
end

if with_RfnXrs && nRefine_inc > nRefine_xrs
    nRefine_xrs = nRefine_inc; % (at least)
end

if with_RfnXrs && ~with_Update
    warning(['adaptive refinement of crack increments requires that ', ...
     '''with_Update'' = TRUE; forcing an override.']); with_Update = 1;
end

%--------------------------------------------------------------------------
% Enrichment parameters
%--------------------------------------------------------------------------

switch mesh_ElemType
    case 'T3'
        
        f_crk = 1.10; % search radius for split elements
        
        switch mesh_EnriSize
            case 'BIG'
                f_tip = 7.00;
                f_sif = 4.00;
            case 'big'
                f_tip = 4.50; % radius for enriching the tip field
                f_sif = 2.50; % radius for evaluating SIF (<f_tip)
            case 'normal'
                f_tip = 2.50;
                f_sif = 1.50;
            case 'small'
                f_tip = 2.00;
                f_sif = 1.50;
            otherwise
                error(['Unknown enrichment size: mesh_EnriSize = ''%s''\n', ...
                    '(n.b. go to error message to see all options)'],mesh_EnriSize)
        end
        
        % for structured meshes
        f_inc = f_tip+1; % (>f_tip) crack increment length
        f_xrs = f_tip+1; % (>f_sif) minimum distance between cracks  
        
        if with_GLwDir
            f_xrs = f_inc+1; % +1 ring of elements
        end
        
        if with_RfnXrs % || with_RfnInc
            f_xrs = f_xrs+1*0.5;
            f_rfn = f_tip + 1.2; % (1.00 < FOS < 2.00)
        end
        
    case 'Q4'
        
        f_crk = 1.10; % search radius for split elements
        
        switch mesh_EnriSize
            case 'HUGE'
                f_tip = 5.99;
                f_sif = f_tip-sqrt(2)*1.50;
            case 'BIG'
                f_tip = 3.99;
                f_sif = f_tip-sqrt(2)*0.50;
            case 'big'
                f_tip = 2.99;               % radius for enriching the tip field
                f_sif = f_tip-sqrt(2)*0.50;  % radius for evaluating SIF (<f_tip) 
            case 'normal'
                f_tip = 2.50; % 1.99;
                f_sif = 1.50;
            case 'small'
                f_tip = 1.50;
                f_sif = f_tip;
            otherwise
                error(['Unknown enrichment size: mesh_EnriSize = ''%s''\n', ...
                    '(n.b. go to error message to see all options)'],mesh_EnriSize)
        end
        
        % for structured meshes
        f_inc = f_tip+sqrt(2); % (>f_tip) crack increment length         
        f_xrs = f_tip+sqrt(2); % (>f_sif) minimum distance between cracks  
        
        if with_GLwDir
            f_xrs = f_inc+sqrt(2); % +1 ring of elements
        end
        
        if with_RfnXrs
            error('Tip refinement is not possible for element type Q4')
        end
        
    otherwise
        error(['Unknown element type: mesh_ElemType = ''%s''\n', ...
            '(n.b. available alternatives are: ''T3'', ''Q4'')'],mesh_ElemType)
end

% param. rel. to f_tip
f_sif = f_sif/f_tip;
f_inc = f_inc/f_tip;
f_xrs = f_xrs/f_tip;
