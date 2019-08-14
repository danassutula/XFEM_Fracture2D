function [ndclr,nddel,eldel] = ...
    Topo_omgStd_remesh(ndcrd_new,elnod_new,elold,ndfre,ndfix,tol)

%--------------------------------------------------------------------------
%
% Topo_omgStd_add:
%   incorporates a re-meshed patch within the existing mesh
%
% INPUT:
%   ndcrd_new   - nodal coordinates of new mesh 
%   elnod_new   - element topology of new mesh 
%   elold       - existing elements that to be replaced with new
%   ndfre       - interior nodes of the existing mesh to be replaced
%   ndfix       - boundary nodes on the existing mesh to be replaced
%   tol         - tolerence for matching nodes on boundary interface
%
% OUTPUT:
%   ndclr - new nodes that have replaced old (to be used to clear eq.)
%   nddel - unused nodes that have been deleted (to be used to delete eq.)
%   eldel - deleted elements (to be used in shifting element related data)
%
%--------------------------------------------------------------------------

global nNdStd mNdCrd nElems mLNodS

nnnew = size(ndcrd_new,1);
nenew = size(elnod_new,1);

neold = length(elold);
nnfre = length(ndfre);
nnfix = length(ndfix);

% find which new nodes correspond to old nodes on boundary
[idmbr,idref] = Nodes_mbr(ndcrd_new,mNdCrd(ndfix,:),tol);

if length(idmbr) ~= nnfix
    
    error('could not find all nodes'); % warning
    ndfre = [ndfre;setdiff(ndfix,ndfix(idref))];
    nnfre = length(ndfre); nnfix = length(idmbr);
    
end

ndclr = ndfre; % init. nd. to clr.
nddel = [];    % init. nd. to del.
eldel = [];    % init. el. to del.

% new nodes to insert/append
idnew = setdiff(1:nnnew,idmbr);

% net n. of nodes added
nnadd = nnnew-(nnfre+nnfix);

% n. of elements added
neadd = nenew-neold;

% backup original topo.
elnod_tmp = elnod_new;

if nnadd > 0
    
    for i = 1:nnfre; k = idnew(i); % reuse
        elnod_new(elnod_tmp==k) = ndfre(i);
        mNdCrd(ndfre(i),:) = ndcrd_new(k,:);
    end
    
    % resize nodal coordiantes
    mNdCrd(end+nnadd,1) = 0;
    
    k = nNdStd; 
    for i = idnew(nnfre+1:end); % append
        k=k+1; elnod_new(elnod_tmp==i) = k;
        mNdCrd(k,:) = ndcrd_new(i,:);
    end
    
else % warning('coarsening mesh (nodes)')
    
    for i = 1:nnfre+nnadd; k = idnew(i);
        elnod_new(elnod_tmp==k) = ndfre(i);
        mNdCrd(ndfre(i),:) = ndcrd_new(k,:);
    end
    
    ndclr = ndfre(1:i);
    nddel = ndfre(i+1:end);
    mNdCrd(nddel,:) = [];
    
end

% update number of nodes
nNdStd = nNdStd + nnadd;

% update number of elements
nElems = nElems + neadd;

% convenient for shifitng topo.
nddel = sort(nddel,'descend');

for i = 1:nnfix % make el. topo. consistent on bnd.
    elnod_new(elnod_tmp==idmbr(i)) = ndfix(idref(i));
end

if neadd > 0
    
    mLNodS(elold,:) = elnod_new(1:neold,:);
    mLNodS(end+1:nElems,:) = elnod_new(neold+1:end,:);
    
else % warning('coarsening mesh (elements)')
    
    mLNodS(elold(1:nenew),:) = elnod_new;
    eldel = elold(nenew+1:end);
    mLNodS(eldel,:) = [];
    
    for i = nddel(:)' % shift element nodes
        q = mLNodS > i; mLNodS(q) = mLNodS(q)-1;
    end

end