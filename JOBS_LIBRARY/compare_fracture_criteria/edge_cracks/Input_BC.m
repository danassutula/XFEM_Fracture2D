

%==========================================================================
% Input Boundary Conditions (BC)
%==========================================================================


%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------

vFxBnd = []; % list of fixed boundaries
mFxDim = []; % bnd. fixed dim.; mFxDim(i,:) = [  1,  0];
mFxVal = []; % bnd. fixed val.; mFxVal(i,:) = [ ux, uy];

vLdBnd = []; % list of loaded boundaries
mBnLod = []; % bnd. traction; mBnLod(i,:) = [f_x,f_y];

mCkLod = []; % internal crack pressure; mCkLod = [p_n,p_t];
vBdLod = []; % body load (volumetric); vBdLod = [b_x,b_y];

hPrLod      = []; % pre-stress/strain (s.t. s = D*(e+e0)+s0) fun. handle
cPrLod_var  = {}; % list (cell) of arguments for the pre-load function

% n.b. concerning pre-load, it is physically meaningful to specify either
% pre-stress (s0) or pre-strain (e0); thus, at least one has to be zero.

% n.b. the listing of domain boundaries, as obtained from the Mesh_make.m
% for a rectangular domain is as follows:
%   1 - bottom edge nodes
%   2 - right edge nodes
%   3 - top edge nodes
%   4 - left edge nodes
%   5 - bottom-left node
%   6 - bottom-right node
%   7 - top-right node
%   8 - top-left node

% (!) The above convention of domain boundary listing can be different if 
% a mesh file is read instead

%--------------------------------------------------------------------------



%==========================================================================
% BC: Dirichlet (essential)
%==========================================================================


%--------------------------------------------------------------------------
% Displacement fixities and values
%--------------------------------------------------------------------------

% vFxBnd = [1,3];
vFxBnd = [5,6];

mFxDim(1,:) = [   1,   1];
mFxVal(1,:) = [-0.0,-0.0];

mFxDim(2,:) = [   0,   1];
mFxVal(2,:) = [ 0.0, 0.0];

%--------------------------------------------------------------------------



%==========================================================================
% BC: Neumann (natural)
%==========================================================================


%--------------------------------------------------------------------------
% Boundary tractions
%--------------------------------------------------------------------------

vLdBnd = [1,3]; % nodes below and above end-crack (RHS)

mBnLod(1,:) = [0,-10];
mBnLod(2,:) = [0, 10];

%--------------------------------------------------------------------------
% Crack surface tractions
%--------------------------------------------------------------------------

wCkLod = 0; % with crack pressure
mCkLod(1,:) = [1,0]; % =[p_n,p_t]

%--------------------------------------------------------------------------
% Body load
%--------------------------------------------------------------------------

wBdLod = 0; % with body load
vBdLod = [0,0]; % =[b_x,b_y]

%--------------------------------------------------------------------------
% Pre-stress/pre-strain/thermal load as a function of position (x,y)
%--------------------------------------------------------------------------

wPrLod = 0; % with pre-load
hPrLod = @Input_BC_preload;
cPrLod_var = {};

%--------------------------------------------------------------------------
