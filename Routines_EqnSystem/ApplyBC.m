function [vFxDof_all,vFxDof_nnz,vFxVal_nnz] = ApplyBC(cFxNod,mFxDim,mFxVal)

% ApplyBC: gets BC DOFs (std)

% INPUT:
%   cFxNod - cell of boundaries that are vectors of fixed nodes, size = {n_fixBound};
%   mFxDim - fixed displacement dimensions of boundaries, size = [n_fixBounds,n_dims], logical (0 or 1);
%   mFxVal - fixed displacement values for boundaries, size = [n_fixBounds,n_dims];

global mNdCrd; [n_nod,n_dim] = size(mNdCrd);
gldof = reshape(1:n_nod*n_dim,n_dim,n_nod);

vFxDof_all = [];
vFxDof_nnz = [];
vFxVal_nnz = [];

for i = 1:length(cFxNod)
    
    fxnod = cFxNod{i};
    fxdim = mFxDim(i,:);
    fxval = mFxVal(i,:).*fxdim;
    
    if size(fxnod,2) ~= 1
        error('Expected a vector of unique nodes for the Dirichlet boundary')
    end
    
    p1 = gldof(fxdim==1,fxnod);
    p2 = gldof(fxval~=0,fxnod);
    p3 = fxval(ones(1,length(fxnod)),fxval~=0)';
    
    vFxDof_all = [ vFxDof_all; p1(:) ];
    vFxDof_nnz = [ vFxDof_nnz; p2(:) ];
    vFxVal_nnz = [ vFxVal_nnz; p3(:) ];
    
end
