
%==========================================================================
% Determine enriched shape functions for enriched elements
%
% N.B.: need to follow the order as in the enr. topology 
%
%==========================================================================

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
cGsEnr_omgShE = cell(nElems,1);
cGsEnr_omgDvE = cell(nElems,1);
%--------------------------------------------------------------------------

% store ck. tip. dir. (useful later)
mTpAlf = zeros(nCrack,2);

for i = 1:nCrack
    x = cCkCrd{i}(1,:)-cCkCrd{i}(2,:);
    mTpAlf(i,1) = atan2(x(2),x(1)); end
for i = 1:nCrack
    x = cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:);
    mTpAlf(i,2) = atan2(x(2),x(1)); end

for i = 1:nCrack
mCkCrd = cCkCrd{i};

%--------------------------------------------------------------------------
% Tip #1 elements
%--------------------------------------------------------------------------
vTpCrd = mCkCrd(1,:);
uTpAlf = mTpAlf(i,1);

% std.
vElBrn = [cBrStd_eTp{i,1};cBrStd_eSp{i,1};cBrStd_eFl{i,1}];

ShapesElm_brnStd;

% bln.
vElBrn = [cBrBln_eSp{i,1};cBrBln_eFl{i,1}];
mBlRmp = [cBrBln_rSp{i,1};cBrBln_rFl{i,1}];

ShapesElm_brnBln;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Tip #2 elements
%--------------------------------------------------------------------------
vTpCrd = mCkCrd(end,:);
uTpAlf = mTpAlf(i,2);

% std.
vElBrn = [cBrStd_eTp{i,2};cBrStd_eSp{i,2};cBrStd_eFl{i,2}];

ShapesElm_brnStd;

% bln.
vElBrn = [cBrBln_eSp{i,2};cBrBln_eFl{i,2}];
mBlRmp = [cBrBln_rSp{i,2};cBrBln_rFl{i,2}];

ShapesElm_brnBln;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Heaviside elements
%--------------------------------------------------------------------------
% std. & bln.
vElHvi = [cHvStd_elm{i};cHvBln_elm{i,1};cHvBln_elm{i,2}];
mHvSgm = [cHvStd_sgm{i};cHvBln_sgm{i,1};cHvBln_sgm{i,2}];

ShapesElm_hviStd;
%--------------------------------------------------------------------------

end