
%--------------------------------------------------------------------------
% Declare Global Variables
%--------------------------------------------------------------------------

% % for cracks
% global cCkCrd
% global mCkJun
% global mTpAct
% global mTpRdi

% % Material
% global cDMatx vElPhz

% Degrees of freedom
global mNDofS mNDspS
global mNDofE mNDspE

% Mesh
global nNdStd mNdCrd
global nElems mLNodS

% Boundary nodes
global nBound cBnNod
global cFxNod cLdNod

% Enrichment tracking
global cLEnDt
global cLNodE
global cXsElm

% Element tracking (Hvi.)
global cHvStd_elm
global cHvStd_sgm 
global cHvBln_elm 
global cHvBln_sgm 
global cHvBln_rmp

% Element tracking (Brn.)
global cBrStd_eTp
global cBrStd_eSp
global cBrStd_eFl 
global cBrBln_eSp
global cBrBln_rSp
global cBrBln_eFl
global cBrBln_rFl

% Enriched node tracking (Hvi.)
global cNdEnr_hvi
global cNdBln_hvi

% Enriched node tracking (Brn.)
global cNdEnr_brn
% global cNdBln_brn % (not used)

% n. of branch enrichment fun.'s
global nEnFun_brn
% global nEnFun_hvi % (not used)

% Gauss point mapping: global, local
global with_MapTyp

% Gauss quadrature (std.)
global nGsStd_omg
global mGsStd_omgShp
global mGsStd_omgDrv
global vGsStd_omgWgt

% Gauss quadrature (enr.)
global cGsEnr_omgShS
global cGsEnr_omgShE
global cGsEnr_omgDvS
global cGsEnr_omgDvE
global cGsEnr_omgWgt

% tolerances
global tol_ref
global tol_abs
