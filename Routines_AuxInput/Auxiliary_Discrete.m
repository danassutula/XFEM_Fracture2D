
%--------------------------------------------------------------------------
% Remove redundant element-paths
%--------------------------------------------------------------------------

% IMPORTANT: use forward slash ('/') in paths
% in order to be platform-safe (e.g. UNIX, PC)

switch mesh_ElemType % safe way of retaining only the used element
    case 'T3'; rmpath(genpath([path_main,'/Routines_Element/Q4']))
    case 'Q4'; rmpath(genpath([path_main,'/Routines_Element/T3']))
    otherwise; error(['Unknown element type: ',mesh_ElemType])
end

%--------------------------------------------------------------------------
% Mesh details
%--------------------------------------------------------------------------

if isempty(mLNodS)
    error('No elements found of specified type: %s',mesh_ElemType)
end

switch mesh_ElemType
    case 'T3' % el. format
        
        mElCrd_loc  = [0,0;1,0;0,1];
        jBNodS      = [1,2;2,3;3,1];
        
    case 'Q4' % el. format
        
        mElCrd_loc  = [-1,-1;1,-1;1,1;-1,1];
        jBNodS      = [1,2;2,3;3,4;4,1];
        
    otherwise
        error(['Unknown element type: ',mesh_ElemType])
end

[nNdStd,nDimes] = size(mNdCrd);
[nElems,nLNodS] = size(mLNodS);
[nElBnd,nBNodS] = size(jBNodS);

nBnFix = length(vFxBnd);
nBnLod = length(vLdBnd);

%--------------------------------------------------------------------------
% Get nodes on domain boundary
%--------------------------------------------------------------------------

vNdBnd = Topo_gamAll(mLNodS,jBNodS); % element edges on domain boundary
cBnNod = Topo_sortBnd(vNdBnd);       % element edges for each boundary
vNdBnd = unique(vNdBnd(:));          % nodes on domain boundary
nBound = length(cBnNod);

%--------------------------------------------------------------------------
% Estimate mesh size
%--------------------------------------------------------------------------

n  = min(100,nElems);    % some samples
j  = randi(nElems,n,1);  % random samples

% init. el. size sample
he_smp = zeros(n,1);

for i = 1:n % random sampling of element size
    he_smp(i) = ElemSize(mNdCrd(mLNodS(j(i),:),:));
end; he_smp = sort(he_smp,'ascend');

% estimate supremum (95% confidence)
he_smp = sum(he_smp)/n+2*std(he_smp);

%--------------------------------------------------------------------------
% Set tip enrichment radii
%--------------------------------------------------------------------------

% enr. el. search dist.
uCkRdi = f_crk*he_smp;

% init. tip enr. radii
mTpRdi = zeros(nCrack,2);
% init. tip activity
mTpAct = false(nCrack,2);
% init. ck. junctions
mCkJun = zeros(nCrack,2);

% est. tip enr. radius
r_tip = f_tip*he_smp;

for i = 1:nCrack
    
    x_tip = cCkCrd{i}(1,:);
    
    q = Elems_tip(mNdCrd,mLNodS,x_tip,r_tip);
    for j = q(:)'; x_elm = mNdCrd(mLNodS(j,:),:);
        if GeoElm_xTip(x_tip,x_elm)
            
            mTpRdi(i,1) = f_tip*ElemSize(x_elm);
            mTpAct(i,1) = true; break
            
        end
    end
    
    x_tip = cCkCrd{i}(end,:);
    
    q = Elems_tip(mNdCrd,mLNodS,x_tip,r_tip);
    for j = q(:)'; x_elm = mNdCrd(mLNodS(j,:),:);
        if GeoElm_xTip(x_tip,x_elm)
            
            mTpRdi(i,2) = f_tip*ElemSize(x_elm);
            mTpAct(i,2) = true; break
            
        end
    end
    
end; clear x_tip x_elm r_tip

% freeze crack tips
for i = 1:size(mTpFrz,1)
    switch mTpFrz(i,2)
        case 1; mTpAct(mTpFrz(i,1),1) = 0;
        case 2; mTpAct(mTpFrz(i,1),2) = 0;
    end
end

% assign same radius for all active tips
if ~with_AdpEnr % && strcmp(kind_GrwCrt,'symmetric')
    i = find(mTpAct); mTpRdi(i) = sum(mTpRdi(i))/length(i);
end

%--------------------------------------------------------------------------
% Set tolerances 
%--------------------------------------------------------------------------

% update reference element size (this is more relevant)
he_ref = mean(reshape(mTpRdi(mTpRdi~=0),[],1))/f_tip;

tol_ref = 1e-12;            % relative to machine (global var.)
tol_abs = tol_ref*he_ref;   % relative to element (global var.)

% N.B.: difficult to determine what the right tolerance should be for 
% elements that are cut by a crack along the edges. This specifically
% regards the following functions: ForceEnr_BndCkHvi, ForceEnr_BndCkBrn, 
% IntJ_CrkLodHvi, IntJ_CrkLodBrn, IntM_CrkLod

%--------------------------------------------------------------------------
% Get BC nodes
%--------------------------------------------------------------------------

cFxNod = []; cFxCrd = []; cFxNod_egs = [];
cLdNod = []; cLdCrd = [];

if nBnFix < size(mFxDim,1) && nBnFix < size(mFxVal,1)
    warning('number of fixed boundaries: n_fix = %d',nBnFix)
    mFxDim = mFxDim(1:nBnFix,:); mFxVal = mFxVal(1:nBnFix,:);
elseif nBnFix > size(mFxDim,1) || nBnFix > size(mFxVal,1)
    error('Undefined fixed boundaries')
end

if nBnLod < size(mBnLod,1)
    warning('number of loaded boundaries: n_lod = %d',nBnLod)
    mBnLod = mBnLod(1:nBnLod,:);
elseif nBnLod > size(mBnLod,1)
    error('Undefined loaded boundaries')
end

if exist('cBCNod','var') && ~isempty(cBCNod)
    
    if max(vFxBnd) > length(cBCNod)
        error('Some fixed boundary(ies) undefined in mesh file');  end
    if max(vLdBnd) > length(cBCNod)
        error('Some loaded boundary(ies) undefined in mesh file'); end
    
    for i = 1:nBnFix
        cFxNod{i,1} = cBCNod{vFxBnd(i)};
    end
    for i = 1:nBnLod
        cLdNod{i,1} = cBCNod{vLdBnd(i)};
    end
    
elseif exist('cBCCrd','var') && ~isempty(cBCCrd)
    
    if max(vFxBnd) > length(cBCCrd)
        error('Some fixed boundary(ies) undefined in mesh file');  end
    if max(vLdBnd) > length(cBCCrd)
        error('Some loaded boundary(ies) undefined in mesh file'); end
    
    for i = 1:nBnFix
        cFxCrd{i,1} = cBCCrd{vFxBnd(i)};
    end
    for i = 1:nBnLod
        cLdCrd{i,1} = cBCCrd{vLdBnd(i)};
    end
    
else
    error('Undefined boundaries')
end

for i = 1:nBnFix
    if isempty(cFxNod{i})
        cFxNod{i,1}=Nodes_mbr(mNdCrd,cFxCrd{i},tol_abs);
        if length(cFxNod{i}) ~= size(cFxCrd{i},1)
            error('Could not find Dirichlet nodes')
        end
    end
end
for i = 1:nBnLod
    if isempty(cLdNod{i})
        cLdNod{i,1}=Nodes_mbr(mNdCrd,cLdCrd{i},tol_abs);
        if length(cLdNod{i}) ~= size(cLdCrd{i},1)
            error('Could not find Neumann nodes')
        end
    end
end

for i = 1:nBnFix
    if size(cFxNod{i},2) == 2
        cFxNod_egs{i} = cFxNod{i}; % [n1,n2;...]
        cFxNod{i} = unique(cFxNod{i}(:)); % [n1,n2,...]'
    else
        if size(cFxNod{i},1) > 1
            cFxNod_egs{i} = Topo_gamSub( ...
                cFxNod{i},mLNodS,jBNodS);
        end
    end
end
for i = 1:nBnLod
    if size(cLdNod{i},2) ~= 2 && size(cLdNod{i},1) > 1 
        % get boundary edge topology from nodes vector 
        cLdNod{i} = Topo_gamSub(cLdNod{i},mLNodS,jBNodS);
    end
end

% no longer required
clear cBCCrd cFxCrd cLdCrd

%--------------------------------------------------------------------------
% Check element phases
%--------------------------------------------------------------------------

if nPhase == 1
    vElPhz = ones(nElems,1);
else
    if ~exist('vElPhz','var') || ...
       length(vElPhz)~=nElems || min(vElPhz)~=1 || max(vElPhz)~=nPhase
        error('Undefined mesh phases; n_phz = %d',nPhase);
    end
end

%--------------------------------------------------------------------------
% Compute element pre-load: S0 = -D*e0 + s0  
%--------------------------------------------------------------------------

if wPrLod % (compute pre-stress/strain at element centers)
    mPrLod = MtxPreload(hPrLod,cPrLod_var,cDMatx,vElPhz,[]);
else
    mPrLod = zeros(3,nElems);
end

%--------------------------------------------------------------------------
% Set not to compute forces if zero load
%--------------------------------------------------------------------------

if wCkLod && ~any(mCkLod(:)); wCkLod = 0; end
if wPrLod && ~any(mPrLod(:)); wPrLod = 0; end
if wBdLod && ~any(vBdLod);    wBdLod = 0; end

%--------------------------------------------------------------------------
% Test atan2 is working as expected
%--------------------------------------------------------------------------

if abs(atan2(0,-1) - pi) > eps
    error(['atan2(Y,X) does not bahave as expected; ',...
     'the code requires that atan2(Y=0,X=-1) == ',num2str(pi),...
     ', but instead the output value was: ',num2str(atan2(0,-1))])
    % matlab implementation of atan2 has changed from the one in version 
    % R2012b the program might not work as expected if cracks are going
    % to be lying on element edges.
end

%--------------------------------------------------------------------------
