function mPrLod = MtxPreload(hPrLod,cPrLod_var,cDMatx,vElPhz,idelm)

%--------------------------------------------------------------------------
%
% mPrLod = pre-load at element centers, size(3,nElems)
%
% hPrLod = pre-load function handle (user defined)
% cPrLod_var = list of function arguments (user defined)
% cDMatx = constitutive matrices for different material phases
% vElPhz = element-wise list of material phases
% idelm  = subset of elements for which pre-load is computed;
%     if 'idelm' is empty - compute for all elements        
%
%--------------------------------------------------------------------------

if nargin ~= 5
    error('Wrong number of input arguments.');
end

global mNdCrd mLNodS; % (CONST)

if isempty(idelm)
    [nelem,nnode]=size(mLNodS);
    idelm = 1:nelem; % all elem.
else
    [nelem,nnode]=size(mLNodS(idelm,:));
    idelm = idelm(:)'; % sub. elem.
end

% eval. pre-load at elm. centers
% full pre-load: S0 = -D*e0 + s0

x_mid = zeros(nelem,2);

for i = 1:nelem
    x_mid(i,:) = sum(mNdCrd(mLNodS(idelm(i),:),:),1);
end

x_mid = x_mid/nnode;

% call pre-load function by handle
[e0,s0] = hPrLod(x_mid,cPrLod_var); 

% build element pre-load matrix
if any(e0(:)~=0) % && ~isempty(e0)
    if min(size(e0))==1; e0 = e0(:);
         mPrLod = e0(:,ones(nelem,1));
    else mPrLod = e0; end
    for i = 1:nelem
        mPrLod(:,i) = cDMatx{vElPhz(idelm(i))}*mPrLod(:,i);
    end
elseif any(s0(:)~=0) % && ~isempty(s0)
    if min(size(s0))==1; s0 = s0(:);
         mPrLod = s0(:,ones(nelem,1));
    else mPrLod = s0; end
else
    mPrLod = zeros(3,nelem);
end

end
