
%==========================================================================
%                                 PART I
%==========================================================================

if ~islogical(vElEnr); p = vElEnr;
    vElEnr = false(nElems,1);
    vElEnr(p) = true;
end

if ~islogical(vElUp0); p = vElUp0;
    vElUp0 = false(nElems,1);
    vElUp0(p) = true;
end

%--------------------------------------------------------------------------
% Elements/nodes in location of mesh refinement
%--------------------------------------------------------------------------

idtip = zeros(nCrack,2); % new crack tip location
q=mCk2Rf(:,1); idtip(q,1)=max(mNoInc_grw(q,1),1);
for i = find(mCk2Rf(:,2))'; n=size(cCkCrd{i},1);
    idtip(i,2) = min(n-mNoInc_grw(i,2)+1,n);
end % (also appropriate for crack junctions)

% tip refinement radii
tprd0 = mTpRdi_ref.*(2.^mCkRfN_ref)*(f_rfn/f_tip);

% init. el. to remove
elems = false(nElems,1);

for i = find(mCk2Rf(:,1))'
    
    elems(Elems_tip(mNdCrd,mLNodS,cCkCrd{i}(idtip(i,1),:), ...
        tprd0(i,1))) = true; % get mesh region at new crack tip
    
    % mesh at prev. tip [to coarsen]
    elems(cElRfn{i,1}) = true;
    
    if ~isempty(cHvBln_elm{i,1}) % Hvi. bln. must be cleared for inc. sgm.
        
        % cast blending elements as standard (Hvi.)
        cHvStd_elm{i} = [cHvBln_elm{i,1};cHvStd_elm{i}];
        cHvStd_sgm{i} = [cHvBln_sgm{i,1};cHvStd_sgm{i}];
        
        cHvBln_elm{i,1} = [];
        cHvBln_sgm{i,1} = [];
        cHvBln_rmp{i,1} = [];
        
    end
    
    cNdBln_hvi{i,1} = [];
    
end

for i = find(mCk2Rf(:,2))'
    
    elems(Elems_tip(mNdCrd,mLNodS,cCkCrd{i}(idtip(i,2),:), ...
        tprd0(i,2))) = true; % get mesh region at new crack tip
 
    % mesh at prev. tip [to coarsen]
    elems(cElRfn{i,2}) = true;
    
    if ~isempty(cHvBln_elm{i,2}) % Hvi. bln. must be cleared for inc. sgm.
        
        % cast blending elements as standard (Hvi.)
        cHvStd_elm{i} = [cHvStd_elm{i};cHvBln_elm{i,2}];
        cHvStd_sgm{i} = [cHvStd_sgm{i};cHvBln_sgm{i,2}];
        
        cHvBln_elm{i,2} = [];
        cHvBln_sgm{i,2} = [];
        cHvBln_rmp{i,2} = [];
        
    end
    
    cNdBln_hvi{i,2} = [];
    
end

while 1 % check if other tips need to be remeshed for mesh comformity
    
    % nodes of current elements to remesh
    % only need boundary nodes; all is faster
    tmp = reshape(mLNodS(elems,:),[],1); % find(elems)
    
    flag = 1; % (break)
    
    for i = setdiff(1:nCrack,find(mCk2Rf(:,1)))
        if any(ismember(tmp,mLNodS(cElRfn{i,1},:))); flag=0;
            elems(cElRfn{i,1}) = true; mCk2Rf(i,1) = true;
            idtip(i,1) = 1; % (okey, not a junction)
        end
    end
    
    for i = setdiff(1:nCrack,find(mCk2Rf(:,2)))
        if any(ismember(tmp,mLNodS(cElRfn{i,2},:))); flag=0;
            elems(cElRfn{i,2}) = true; mCk2Rf(i,2) = true;
            idtip(i,2) = size(cCkCrd{i},1); % (also okey)
        end
    end
    
    if flag
        break % no other tips to remesh
    end
    
end

elems = find(elems);
nelem = length(elems);

vTp2Rf = find(mCk2Rf(:));
nTp2Rf = length(vTp2Rf);

% refine tips in order of coarsest to finest
[~,q] = sort(mCkRfN(vTp2Rf)); vTp2Rf=vTp2Rf(q);
[q_crk,q_tip] = ind2sub([nCrack,2],vTp2Rf);

%--------------------------------------------------------------------------
% Elements in location of mesh refinement
%--------------------------------------------------------------------------

q_up0 = ismember(elems,find(vElUp0)); elup0 = elems(q_up0); % enr. el. whose enr. eq. were rmv.
q_std =~ismember(elems,find(vElEnr)); elstd = elems(q_std); % std. el. that have std. eq. only
q_upd =~(q_std|q_up0);                elupd = elems(q_upd); % enr. el. that have all their eq.

%--------------------------------------------------------------------------
% Nodes on the interface with the rest of the mesh 
%--------------------------------------------------------------------------

ndsd0 = Topo_gamAll(mLNodS(elems,:),jBNodS);
ndsd0 = unique(ndsd0(:)); % interface nodes

%--------------------------------------------------------------------------
% Discard nodes that happen to lie on domain boundary
%--------------------------------------------------------------------------

ndbnd = intersect(ndsd0,vNdBnd);
if ~isempty(ndbnd) % lie on bnd.
    
    % get element edges that contain nodes ndbnd
    tmp = Topo_gamSub(ndbnd,mLNodS(elems,:),jBNodS);
    
    % (el. lying on ndbnd will be remeshed)
    ndbnd = intersect(tmp(:,1),tmp(:,2));
    ndsd0 = setdiff(ndsd0,ndbnd); % keep
    
end

%--------------------------------------------------------------------------
% Nodes in the interior of mesh refinement
%--------------------------------------------------------------------------

ndstd = unique(reshape(mLNodS(elems,:),[],1));
ndstd = setdiff(ndstd,ndsd0); % interior nodes
% nddel = []; % (get later)

nnsd0 = length(ndsd0);
nnstd = length(ndstd);

%--------------------------------------------------------------------------
% Disassemble global stiffness
%   - sufficient to decouple a patch from global system of equations
%   - all interior degrees of freedom (std/enr) will be deleted later
%--------------------------------------------------------------------------

% elup0 (\in vElUp0) elements have no enriched shapes, only the leftover
% enriched nodes; the enriched contribution to Kg, Fg already subtracted

% enough to rmv. these el. so interior is decoupled
elbln = elems(any(ismember(mLNodS(elems,:),ndsd0),2));

% group as std. since enr. in elup0 already removed
elstd = [elstd;elup0];

for i = intersect(elstd,elbln)' % (decoupling)
    
    q = nDimes*mLNodS(i,:);
    q = [q-1;q]; q=q(:);
    
    % subtract Ke
    mGlStf(q,q) = mGlStf(q,q) - MtxStiffElm(mNdCrd(mLNodS(i,:),:),...
        cDMatx{vElPhz(i)},mGsStd_omgDrv,vGsStd_omgWgt);
    
    % subtract Fe (body)
    if wBdLod
        vGlFrc(q) = vGlFrc(q) - ForceElm_body(mNdCrd(mLNodS(i,:),:),...
            mGsStd_omgShp,mGsStd_omgDrv,vGsStd_omgWgt,vBdLod);
    end
    
    % subtract Fe (pre-load)
    if wPrLod
        vGlFrc(q) = vGlFrc(q) - ForceElm_resid(mNdCrd(mLNodS(i,:),:),...
            mGsStd_omgDrv,vGsStd_omgWgt,mPrLod(:,i));
    end
    
end

% elupd (\in vElEnr) elements have enriched shapes and nodes;
% the enriched and standard contributions to Kg, Fg are subtracted

for i = intersect(elupd,elbln)' % (decoupling)
    
    q = nDimes*mLNodS(i,:);
    q = [q-1;q];
    
    p = mNDofE(:,cLNodE{i});
    q = [q(:);p(:)];
    
    % subtract Ke
    mGlStf(q,q) = mGlStf(q,q) - MtxStiffElm(mNdCrd(mLNodS(i,:),:),...
        cDMatx{vElPhz(i)},[cGsEnr_omgDvS{i},cGsEnr_omgDvE{i}],cGsEnr_omgWgt{i});
    
    % subtract Fe (body)
    if wBdLod
        vGlFrc(q) = vGlFrc(q) - ForceElm_body(mNdCrd(mLNodS(i,:),:),...
            [cGsEnr_omgShS{i},cGsEnr_omgShE{i}],cGsEnr_omgDvS{i},...
            cGsEnr_omgWgt{i},vBdLod);
    end
    
    % subtract Fe (pre-load)
    if wPrLod
        vGlFrc(q) = vGlFrc(q) - ForceElm_resid(mNdCrd(mLNodS(i,:),:),...
            [cGsEnr_omgDvS{i},cGsEnr_omgDvE{i}],cGsEnr_omgWgt{i},mPrLod(:,i));
    end
    
end

%--------------------------------------------------------------------------
% Disasseble enriched topology
%--------------------------------------------------------------------------

% all enr. el. in zone
elenr = [elup0;elupd];

% does everything (n.b. ndenr is logical)
[ndenr,enunq] = Topo_omgEnr_rmv(elenr);
if isempty(enunq); enunq=zeros(0,2); end

% assert Hvi. enr. for grown cracks is included
tmp = find(any(mCk2Up & mCkRfn,2)); tmp(end,2)=0;

enunq = union(enunq,tmp,'rows');
q = enunq(:,2)==0; % branch first
enunq = [enunq(~q,:);enunq(q,:)];

for i = elupd' % = elenr XOR elup0
    cGsEnr_omgShS{i} = [];
    cGsEnr_omgDvS{i} = [];
    cGsEnr_omgDvE{i} = [];
    cGsEnr_omgShE{i} = [];
    cGsEnr_omgWgt{i} = [];
end

% nolonger enriched
vElEnr(elenr) = false;
vElUp0(elup0) = false;

%--------------------------------------------------------------------------
% Refine patch
%--------------------------------------------------------------------------

% initialize replacement el.
elem0 = false(nElems_ref,1);

% backup new and old el. for the relevant tips
for i = 1:nTp2Rf; i_crk=q_crk(i); i_tip=q_tip(i);
    
    if i_tip==1; x_tip=cCkCrd{i_crk}(idtip(i_crk,1),:);
    else x_tip=cCkCrd{i_crk}(idtip(i_crk,2),:); end
    
    % el. from prev. tip refinement
    elem0(cElRf0{i_crk,i_tip}) = 1;
    
    if mCk2Up(i_crk,i_tip) % get el. at current tip
        elem0(Elems_tip(mNdCrd_ref,mLNodS_ref, ...
            x_tip,tprd0(i_crk,i_tip))) = 1;
    end
    
end

elem0 = find(elem0);

elnod_new = mLNodS_ref(elem0,:); % el. topology to remesh
ndcrd_new = mNdCrd_ref(unique(elnod_new(:)),:); % (ascending)
[~,elnod_new] = Topo_omgRef(0,elnod_new); % (reference topo.)

% id of tips with refin. from coarsest to finest
rfunq = unique(mCkRfN(vTp2Rf),'first'); k = 1;

for i = 1:length(rfunq) % (sequential remeshing)
    
    % tips with equal or greater refinement level
    q = vTp2Rf(mCkRfN(vTp2Rf) >= rfunq(i)); q=q(:)';
    
    while k <= rfunq(i) % keep refining
        
        % el. to refine (reg. split - not bisection)
        p = false(size(elnod_new,1),1);
        
        for j = q
            if j > nCrack; j=j-nCrack;
                p(Elems_tip(ndcrd_new,elnod_new, ...
                 cCkCrd{j}(idtip(j,2),:),tprd0(j,2)/2^(k-1))) = true;
            else
                p(Elems_tip(ndcrd_new,elnod_new, ...
                 cCkCrd{j}(idtip(j,1),:),tprd0(j,1)/2^(k-1))) = true;
            end
        end
        
        if ~any(p)
            error('No elements could be found for splitting')
        end
        
        if k == 1
            if isempty(ndbnd)
                
                % find elements to be split; these are the elements that
                % contain no nodes on the patch boundary; refinement of
                % elements on the patch boundary must be prohibited
                
                p(any(ismember(elnod_new,unique(Topo_gamAll(...
                    elnod_new(p,:),jBNodS))),2)) = false;
                
            else
                
                % find elements to be split; these are the elements that
                % contain no nodes on the patch boundary excluding those
                % nodes on the patch boundary that coincide with the nodes
                % of the domain boundary; refinement is allowed on domain 
                % boundary
                %
                % "unique(Topo_gamAll(elnod_new(p,:),jBNodS))" are all the
                % nodes on the boundary of the new patch to be refined;
                %
                % Nodes_mbr(ndcrd_new,mNdCrd(ndbnd,:),tol_abs)) are the
                % corresponding member nodes of the new patch; these nodes
                % lie on the boundary of the domain
                
                p(any(ismember(elnod_new,setdiff(unique(Topo_gamAll(elnod_new(p,:),jBNodS)),...
                    Nodes_mbr(ndcrd_new,mNdCrd(ndbnd,:),tol_abs))),2)) = false;
                
            end
        end
        
        % simple element division: bisection and regular (1 to 4)
        [ndcrd_new,elnod_new] = refine(ndcrd_new,elnod_new,p);
        
        % el. to smooth (safer selection)
        p = false(size(elnod_new,1),1);
        
        for j = q
            if j > nCrack; j=j-nCrack;
                p(Elems_tip(ndcrd_new,elnod_new, ...
                 cCkCrd{j}(idtip(j,2),:),tprd0(j,2))) = true;
            else 
                p(Elems_tip(ndcrd_new,elnod_new, ...
                 cCkCrd{j}(idtip(j,1),:),tprd0(j,1))) = true;
            end
        end
        
        % Laplacian mesh smoothing (boundary nodes are kept fixed)
        ndcrd_new = MeshSmooth(ndcrd_new,elnod_new(p,:),jBNodS);
        
        k = k + 1;
        
    end
end

%--------------------------------------------------------------------------
% Update mesh
%--------------------------------------------------------------------------

% nnode_new = size(ndcrd_new,1);
% nelem_new = size(elnod_new,1);

nElem0 = nElems; % backup
nNdSd0 = nNdStd; % backup

% update mesh
[ndclr,nddel,eldel] = Topo_omgStd_remesh(...
 ndcrd_new,elnod_new,elems,ndstd,ndsd0,tol_abs);

% (nddel is already sorted in descending order)
% if ~issorted(nddel(end:-1:1))
%     nddel = sort(nddel,'descend');
% end

% (all updated implicitly)
% nNdStd = size(mNdCrd,1);
% nElems = size(mLNodS,1); 

ndclr = ndclr(:);
nddel = nddel(:);
eldel = eldel(:);

nnadd = nNdStd-nNdSd0; % nd. added
nedel = length(eldel); % el. deleted 

% get new el. (based on how new topology was obtained)
vElNew = [elems(1:end-nedel);(nElem0-nedel+1:nElems)'];

% backup current and reference el. for remeshed tips
for i = 1:nTp2Rf; i_crk=q_crk(i); i_tip=q_tip(i);
    if mCkRfN(i_crk,i_tip) % backup only if needed
        
        if i_tip==1; x_tip = cCkCrd{i_crk}(idtip(i_crk,1),:);
        else x_tip = cCkCrd{i_crk}(idtip(i_crk,2),:); end
        
        cElRfn{i_crk,i_tip} = vElNew(Elems_tip( ... % current elements
            mNdCrd,mLNodS(vElNew,:),x_tip,tprd0(i_crk,i_tip)));
        
        cElRf0{i_crk,i_tip} = elem0(Elems_tip( ... % reference elements
            mNdCrd_ref,mLNodS_ref(elem0,:),x_tip,tprd0(i_crk,i_tip)));
        
    else % clear redundant elements
        cElRfn{i_crk,i_tip} = [];
        cElRf0{i_crk,i_tip} = [];
    end
end

if nedel > 0 % remove element gaps
   
    q = setdiff(1:nElem0,eldel);
    
    vElEnr = vElEnr(q); nElEnr = nnz(vElEnr);
    vElUp0 = vElUp0(q); nElUp0 = nnz(vElUp0);
    
    cLNodE = cLNodE(q);
    cLEnDt = cLEnDt(q);
    cXsElm = cXsElm(q,:);
    
    cGsEnr_omgShS = cGsEnr_omgShS(q);
    cGsEnr_omgDvS = cGsEnr_omgDvS(q);
    cGsEnr_omgShE = cGsEnr_omgShE(q);
    cGsEnr_omgDvE = cGsEnr_omgDvE(q);
    cGsEnr_omgWgt = cGsEnr_omgWgt(q);
    
elseif nElems > nElem0 % resize for new el.
    
    vElEnr(nElems) = false;
    vElUp0(nElems) = false;
    
    cLNodE{nElems,1} = [];
    cLEnDt{nElems,1} = [];
    cXsElm{nElems,2} = [];
    
    cGsEnr_omgShS{nElems,1} = [];
    cGsEnr_omgDvS{nElems,1} = [];
    cGsEnr_omgShE{nElems,1} = [];
    cGsEnr_omgDvE{nElems,1} = [];
    cGsEnr_omgWgt{nElems,1} = [];
    
end

if ~isempty(eldel) % shift elements if some were deleted
    
    % get previously refined patches
    q = setdiff(find(mCkRfn(:)),vTp2Rf);
    
    if ~isempty(q)
        for i = q(:)'
            for j = 1:length(cElRfn{i}); cElRfn{i}(j)=...
                cElRfn{i}(j) - nnz(cElRfn{i}(j) > eldel);
            end
        end
    end
    
    for i = 1:nCrack
        % Branch el.
        for j = 1:length(cBrStd_eTp{i,1}); cBrStd_eTp{i,1}(j) = ...
                cBrStd_eTp{i,1}(j) - nnz(cBrStd_eTp{i,1}(j) > eldel); end
        for j = 1:length(cBrStd_eTp{i,2}); cBrStd_eTp{i,2}(j) = ...
                cBrStd_eTp{i,2}(j) - nnz(cBrStd_eTp{i,2}(j) > eldel); end
        for j = 1:length(cBrStd_eSp{i,1}); cBrStd_eSp{i,1}(j) = ...
                cBrStd_eSp{i,1}(j) - nnz(cBrStd_eSp{i,1}(j) > eldel); end
        for j = 1:length(cBrStd_eSp{i,2}); cBrStd_eSp{i,2}(j) = ...
                cBrStd_eSp{i,2}(j) - nnz(cBrStd_eSp{i,2}(j) > eldel); end
        for j = 1:length(cBrStd_eFl{i,1}); cBrStd_eFl{i,1}(j) = ...
                cBrStd_eFl{i,1}(j) - nnz(cBrStd_eFl{i,1}(j) > eldel); end
        for j = 1:length(cBrStd_eFl{i,2}); cBrStd_eFl{i,2}(j) = ...
                cBrStd_eFl{i,2}(j) - nnz(cBrStd_eFl{i,2}(j) > eldel); end
        for j = 1:length(cBrBln_eSp{i,1}); cBrBln_eSp{i,1}(j) = ...
                cBrBln_eSp{i,1}(j) - nnz(cBrBln_eSp{i,1}(j) > eldel); end
        for j = 1:length(cBrBln_eSp{i,2}); cBrBln_eSp{i,2}(j) = ...
                cBrBln_eSp{i,2}(j) - nnz(cBrBln_eSp{i,2}(j) > eldel); end
        for j = 1:length(cBrBln_eFl{i,1}); cBrBln_eFl{i,1}(j) = ...
                cBrBln_eFl{i,1}(j) - nnz(cBrBln_eFl{i,1}(j) > eldel); end
        for j = 1:length(cBrBln_eFl{i,2}); cBrBln_eFl{i,2}(j) = ...
                cBrBln_eFl{i,2}(j) - nnz(cBrBln_eFl{i,2}(j) > eldel); end
        % Heaviside el.
        for j = 1:length(cHvStd_elm{i});   cHvStd_elm{i}(j)   = ...
                cHvStd_elm{i}(j)   - nnz(cHvStd_elm{i}(j)   > eldel); end
        for j = 1:length(cHvBln_elm{i,1}); cHvBln_elm{i,1}(j) = ...
                cHvBln_elm{i,1}(j) - nnz(cHvBln_elm{i,1}(j) > eldel); end
        for j = 1:length(cHvBln_elm{i,2}); cHvBln_elm{i,2}(j) = ...
                cHvBln_elm{i,2}(j) - nnz(cHvBln_elm{i,2}(j) > eldel); end
    end
    
end

if ~isempty(ndbnd) % redo domain boundary
    
    vNdBnd = Topo_gamAll(mLNodS,jBNodS); % (edge nodes, N-by-2)
    cBnNod = Topo_sortBnd(vNdBnd);
    vNdBnd = unique(vNdBnd(:));
    
    for i = sort(ndbnd,'descend')'
        for j = 1:nBnFix
            if any(cFxNod{j} == i)
                warning('Dirichlet boundary was remeshed; need to redo BC.');
                % REDO BOUNDARY (cBCNod,cFxNod,cLdNod,cBnNod,vNdBnd):
                % 1) need to use a level set to determine which new nodes
                % are on the old boundary, so that the specific boundary
                % conditions can be prescribed
                % 2) nddel then should not be used to do any shifting of
                % boundary nodes (see below to understand what is ment)  
            end
        end
        for j = 1:nBnLod
            if any(cLdNod{j} == i)
                warning('Neumann boundary was remeshed; need to redo BC.');
                % REDO BOUNDARY:
                % ...
            end
        end
    end
else
    % shift BC nodes if mesh is coarsined
    for i = nddel' % sort(nddel,'descend')'
        
        p=vNdBnd>i; vNdBnd(p)=vNdBnd(p)-1;
        for j = 1:nBound; p = cBnNod{j}>i;
            cBnNod{j}(p) = cBnNod{j}(p)-1;
        end
        
        % alt. for updating vNdBnd:
        % vNdBnd = cat(1,cBnNod{:});
        % vNdBnd = vNdBnd(:,1);
        
    end
end

% IF BC-BOUNDARY IS KEPT THE SAME!

% shift BC nodes if mesh is coarsined
for i = nddel' % sort(nddel,'descend')'
    for j = 1:nBnFix; p = cFxNod{j}>i;
        cFxNod{j}(p) = cFxNod{j}(p)-1;
    end
    
    for j = 1:nBnLod; p = cLdNod{j}>i;
        cLdNod{j}(p) = cLdNod{j}(p)-1;
    end
    
    % shift fixed nodes/dofs
    p = vFxDof_all > nDimes*i;
    vFxDof_all(p) = vFxDof_all(p)-nDimes;
    
    % shift fixed nodes/dofs (non-zero)
    p = vFxDof_nnz > nDimes*i;
    vFxDof_nnz(p) = vFxDof_nnz(p)-nDimes;
end

%--------------------------------------------------------------------------
% Renew element material phases
%--------------------------------------------------------------------------

if nPhase == 1
    if nElems >= nElem0
        vElPhz(end+1:nElems) = vElPhz(1);
    else
        vElPhz = vElPhz(1:nElems);
    end
else % nPhase > 1
    
    % phazes that are affected
    tmp = unique(vElPhz(elems));
    
    if length(tmp) == 1
        if nElems >= nElem0
            vElPhz(end+1:nElems) = tmp;
        else
            vElPhz(eldel) = [];
        end
    else
        
        fStop = 1; % execution will be stopped at the end of 'iStep'
        fStop_log.dbg = dbstack('-completenames'); % debug information
        fStop_log.msg = sprintf('Cannot remesh in multiple materials.');
        
        warning(fStop_log.msg)
        
        % when remeshing takes place in multiple material phases, it is
        % more difficult to determine which (new) elements belong to which
        % material phase; in these circumstances it is necessary (probably)
        % to use the level set approach to sort the elements appropriately
        % or, more simply, derive a function, similar to 'Input_BC_preload'
        % , to determine material phases based on switch/if statements
        
        q = tmp;
        
        for i = 1:length(q) % get dominant material phase
            q(i) = sum(ismember(vElPhz(elems),q(i)));
        end
        
        [~,i]=max(q); tmp=tmp(i);
        vElPhz(elems) = tmp;
        
        if nElems >= nElem0
            vElPhz(end+1:nElems) = tmp;
        else
            vElPhz(eldel) = [];
        end
        
    end
end

%--------------------------------------------------------------------------
% Renew element pre-load
%--------------------------------------------------------------------------

% need to resize even if zeros
if nElems >= nElem0
    mPrLod(3,nElems) = 0;
else
    mPrLod(:,eldel) = [];
end

if wPrLod % (compute pre-stress/strain at element centers)
    mPrLod(:,vElNew) = MtxPreload(hPrLod,cPrLod_var,cDMatx,vElPhz,vElNew);
end

%--------------------------------------------------------------------------
% Clear/delete DOFs
%--------------------------------------------------------------------------

% del. enr. dofs before std. dofs (because: ndenr > nddel)
 
% get dofs to del.
q = mNDofE(:,ndenr);
q = q(:);

% del. dofs from eq.
mGlStf(:,q) = [];
mGlStf(q,:) = [];
vGlFrc(q,:) = [];

% clr. std. dofs before del. std. dofs
 
% get dofs to clr.
q = nDimes*ndclr; 
q = [q-1;q];
 
% clr. dofs from eq.
mGlStf(:,q) = 0;
mGlStf(q,:) = 0;
vGlFrc(q,:) = 0;

% okey, now del. std. dofs

q = nDimes*nddel; 
q = [q-1;q];

% del. dofs from eq.
mGlStf(:,q) = [];
mGlStf(q,:) = [];
vGlFrc(q,:) = [];

%--------------------------------------------------------------------------
% Reclassify enriched DOFs and resize the system of equations
%--------------------------------------------------------------------------

% null cleared dofs
mNDofE(:,ndenr) = 0;

% update enriched nodes
vNdEnr = find(mNDofE(1,:));
nNdEnr = length(vNdEnr);
nGDofE = nDimes*nNdEnr;

nGDofS = nDimes*nNdStd;
nGlDof = nGDofS+nGDofE;

for i = 1:nDimes % shift dofs since some were removed (tmp)
    mNDofE(i,vNdEnr) = nGDofS+i:nDimes:nGlDof-nDimes+i;
end

if nnadd > 0 % resize eq. (make space for any new nodes)
    
    mGlStf(nGlDof,nGlDof) = 0;
    vGlFrc(nGlDof,1)      = 0;
    
    % push down enr. dofs/eq.'s
    q = nGDofS+1:nGlDof;
    n = nDimes*nnadd;
    p = q-n;
    
    % shift Kes, Kse, Kee, Fe
    mGlStf(:,q) = mGlStf(:,p);
    mGlStf(q,:) = mGlStf(p,:);
    vGlFrc(q,:) = vGlFrc(p,:);
    
    % null eq.'s for new std. dofs
    q = nGDofS-n+1:nGDofS;
    
    % zero old values
    mGlStf(:,q) = 0;
    mGlStf(q,:) = 0;
    vGlFrc(q,:) = 0;

end

%--------------------------------------------------------------------------



%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Safety checks
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if with_Debug
    
    elbad_ch1 = [];
    elbad_ch2 = [];
    elbad_ch3 = [];
    
    for i = find(vElEnr(:))'
        if any(reshape(~mNDofE(:,cLNodE{i}),[],1))
            elbad_ch1 = [elbad_ch1;i];
        end
        if isempty(cLNodE{i}) || isempty(cLEnDt{i})
            elbad_ch2 = [elbad_ch1;i];
        end
        if size(unique(cLEnDt{i}(1:2,:)','rows'),1) ~= size(cLEnDt{i},2)
            elbad_ch3 = [elbad_ch3;i];
        end
    end
    
    fail = 0;
    
    if nGlDof ~= length(vGlFrc)
        warning('ch1')
        fail = 1;
    end
    if ~isempty(elbad_ch1)
        warning('ch2')
        fail = 1;
    end
    if ~isempty(elbad_ch2)
        warning('ch3')
        fail = 1;
    end
    if ~isempty(elbad_ch3)
        warning('ch4')
        fail = 1;
    end
    if isempty(enunq) % no enrichments in the remeshing zone
        warning('ch5')
        fail = 1;
    end
    
    if fail
        error('safety check(s) failed')
    end

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%--------------------------------------------------------------------------
% Assemble global stiffness matrix (std)
%--------------------------------------------------------------------------

mGlStf = mGlStf + MtxStiffStd(nGlDof,mNdCrd,mLNodS(vElNew,:),...
    vElPhz(vElNew),cDMatx,mGsStd_omgDrv,vGsStd_omgWgt);

%--------------------------------------------------------------------------
% Boundary force (std)
%--------------------------------------------------------------------------

% ...

%--------------------------------------------------------------------------
% Body tractions (std)
%--------------------------------------------------------------------------
if wBdLod
    
    vGlFrc = vGlFrc + ForceStd_Body(nGlDof,mNdCrd,mLNodS(vElNew,:),...
        mGsStd_omgShp,mGsStd_omgDrv,vGsStd_omgWgt,vBdLod);
    
end
%--------------------------------------------------------------------------
% Pre-stress/pre-strain/thermal (std)
%--------------------------------------------------------------------------
if wPrLod
    
    vGlFrc = vGlFrc + ForceStd_Resid(nGlDof,mNdCrd,mLNodS(vElNew,:),...
        mGsStd_omgDrv,vGsStd_omgWgt,mPrLod(:,vElNew));
    
end
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Delete "local" variables
%--------------------------------------------------------------------------

clear idtip tprd0 q_tp1 q_tp2 q_crk q_tip q_up0 q_std q_upd elems elem0 elbln ...
 elstd elenr eldel ndsd0 ndstd ndenr ndbnd ndclr nddel rfunq ndcrd_new elnod_new

%--------------------------------------------------------------------------




%==========================================================================
%                                 PART II
%==========================================================================




%--------------------------------------------------------------------------
% Add old enrichment
%--------------------------------------------------------------------------

for i_enr = 1:size(enunq,1)
    
    i_crk = enunq(i_enr,1);
    mCkCrd = cCkCrd{i_crk};
    
    switch enunq(i_enr,2)
        case 1 % (brn. tip #1)
        
        %------------------------------------------------------------------
        % Get elements
        % 
        %   - only concerned about this tip enrichment
        %   - ignore other enrichments associated with this crack
        %   - adjust enrichment radius based on refined mesh
        %
        %------------------------------------------------------------------
        
        % need for completeness because of the way elements are determined
        p = unique([vElNew;cBrStd_eTp{i_crk,1};cBrStd_eSp{i_crk,1};...
         cBrStd_eFl{i_crk,1};cBrBln_eSp{i_crk,1};cBrBln_eFl{i_crk,1}]);   
        
        [~,~,~,~,~,eltip_std,elspl_std,elful_std,elspl_bln,elful_bln,elspl_rmp,elful_rmp] ...
         = Elems_enr(mNdCrd,mLNodS(p,:),mCkCrd(1:2,:),[0,0],[mTpRdi(i_crk,1),0],he_ref);
        
        if isempty(eltip_std) % warning
            
            mTpRdi(i_crk,1) = 0;
            mTpAct(i_crk,1) = 0;
            
            error('no tip element?')
            break
        
        end
        
        eltip_std = p(eltip_std{1}); cBrStd_eTp{i_crk,1} = eltip_std;
        elspl_std = p(elspl_std{1}); cBrStd_eSp{i_crk,1} = elspl_std;
        elful_std = p(elful_std{1}); cBrStd_eFl{i_crk,1} = elful_std;
        elspl_bln = p(elspl_bln{1}); cBrBln_eSp{i_crk,1} = elspl_bln;
        elful_bln = p(elful_bln{1}); cBrBln_eFl{i_crk,1} = elful_bln;
        
        cBrBln_rSp{i_crk,1} = elspl_rmp{1};
        cBrBln_rFl{i_crk,1} = elful_rmp{1};
        
        tmp = [eltip_std;elspl_std;elful_std;elspl_bln;elful_bln];
        
        vElBrn = vElNew(ismember(vElNew,tmp));
        vElBr0 = tmp(~ismember(tmp,vElBrn));
        
        vElEnr(vElBrn) = true; % el. is enr.
        vElUp0(vElBrn) = true; % el. to upd.
        
        %------------------------------------------------------------------
        % Update enr. el. topology (glue new to old)
        %------------------------------------------------------------------
        
        [nNdEn8,mLNodE] = ... % get enr. topology (disjoint)
            Topo_omgEnr(nNdEn8,nEnFun_brn,mLNodS(vElBrn,:));
        
        [nNdEn8,mLNodE] = ... % "glue" new topo. to existing
            Topo_omgEnr_add(nNdEn8,mLNodE,vElBrn,vElBr0,i_crk,1);
        
        for i = 1:length(vElBrn) % store enr. nodes (push back)
            cLNodE{vElBrn(i)}(end+1:end+nLNodE_brn) = mLNodE(i,:);
        end
        
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        
        cNdEnr_brn{i_crk,1} = ...
            [cNdEnr_brn{i_crk,1};(nNdEn0+1:nNdEn8)']; nNdEn0 = nNdEn8;
        
        %------------------------------------------------------------------
        % Save el. enr. data: ck. id, enr. id, n. enr. node, GP.
        %------------------------------------------------------------------
        q = ismember(vElBrn,eltip_std);
        
        for i = vElBrn( q)' % tip element(s)
            cLEnDt{i}(:,end+1) = [ i_crk ; 1 ; nLNodE_brn ; 3]; end
        for i = vElBrn(~q)' % other element(s)
            cLEnDt{i}(:,end+1) = [ i_crk ; 1 ; nLNodE_brn ; 2]; end
        
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        
        x_sgm = mCkCrd([1,2],:); % set tip first
        
        for i = eltip_std'
            
            [n,x] = GeoXElem(x_sgm,mNdCrd(mLNodS(i,:),:));
            
            if n == 2
                cXsElm{i,2}(end+1,:) = [1,2]+size(cXsElm{i,1},1);
            end
            
            cXsElm{i,1}(end+1:end+n,:) = x;
            
        end
        
        %------------------------------------------------------------------
        
        case 2 % (brn. tip #2)
        
        %------------------------------------------------------------------
        % Get elements
        %------------------------------------------------------------------
        
        % need for completeness because of the way elements are determined
        p = unique([vElNew;cBrStd_eTp{i_crk,2};cBrStd_eSp{i_crk,2};...
         cBrStd_eFl{i_crk,2};cBrBln_eSp{i_crk,2};cBrBln_eFl{i_crk,2}]);    
        
        [~,~,~,~,~,eltip_std,elspl_std,elful_std,elspl_bln,elful_bln,elspl_rmp,elful_rmp] ...
         = Elems_enr(mNdCrd,mLNodS(p,:),mCkCrd(end-1:end,:),[0,0],[0,mTpRdi(i_crk,2)],he_ref);
        
        if isempty(eltip_std) % warning
            
            mTpRdi(i_crk,2) = 0;
            mTpAct(i_crk,2) = 0;
            
            error('no tip element?')
            break
            
        end
        
        eltip_std = p(eltip_std{2}); cBrStd_eTp{i_crk,2} = eltip_std;
        elspl_std = p(elspl_std{2}); cBrStd_eSp{i_crk,2} = elspl_std;
        elful_std = p(elful_std{2}); cBrStd_eFl{i_crk,2} = elful_std;
        elspl_bln = p(elspl_bln{2}); cBrBln_eSp{i_crk,2} = elspl_bln;
        elful_bln = p(elful_bln{2}); cBrBln_eFl{i_crk,2} = elful_bln;
        
        cBrBln_rSp{i_crk,2} = elspl_rmp{2};
        cBrBln_rFl{i_crk,2} = elful_rmp{2};
        
        tmp = [eltip_std;elspl_std;elful_std;elspl_bln;elful_bln];
        
        vElBrn = vElNew(ismember(vElNew,tmp));
        vElBr0 = tmp(~ismember(tmp,vElBrn));
        
        vElEnr(vElBrn) = true; % el. is enr.
        vElUp0(vElBrn) = true; % el. to upd.
        
        %------------------------------------------------------------------
        % Update enr. el. topology (glue new to old)
        %------------------------------------------------------------------
        
        [nNdEn8,mLNodE] = ... % get enr. topology (disjoint)
            Topo_omgEnr(nNdEn8,nEnFun_brn,mLNodS(vElBrn,:));
        
        [nNdEn8,mLNodE] = ... % "glue" new topo. to existing
            Topo_omgEnr_add(nNdEn8,mLNodE,vElBrn,vElBr0,i_crk,2);
        
        for i = 1:length(vElBrn) % store enr. nodes (push back)
            cLNodE{vElBrn(i)}(end+1:end+nLNodE_brn) = mLNodE(i,:);
        end
        
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        
        cNdEnr_brn{i_crk,2} = ...
            [cNdEnr_brn{i_crk,2};(nNdEn0+1:nNdEn8)']; nNdEn0 = nNdEn8;
        
        %------------------------------------------------------------------
        % Save el. enr. data: ck. id, enr. id, n. enr. node, GP.
        %------------------------------------------------------------------
        q = ismember(vElBrn,eltip_std);
        
        for i = vElBrn( q)' % tip element(s)
            cLEnDt{i}(:,end+1) = [ i_crk ; 2 ; nLNodE_brn ; 3]; end
        for i = vElBrn(~q)' % other element(s)
            cLEnDt{i}(:,end+1) = [ i_crk ; 2 ; nLNodE_brn ; 2]; end
        
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        
        x_sgm = mCkCrd([end,end-1],:); % set tip first
        
        for i = eltip_std'
            
            [n,x] = GeoXElem(x_sgm,mNdCrd(mLNodS(i,:),:));
            
            if n == 2
                cXsElm{i,2}(end+1,:) = [1,2]+size(cXsElm{i,1},1);
            end
            
            cXsElm{i}(end+1:end+n,:) = x;
            
        end
        
        %------------------------------------------------------------------
        
        otherwise % == 0 (Hvi.)
        
        %------------------------------------------------------------------
        % Get elements
        %------------------------------------------------------------------
        
        % remaining Heaviside-enriched elements for this crack 
        vElHv0=[cHvBln_elm{i_crk,1};cHvStd_elm{i_crk};cHvBln_elm{i_crk,2}];
        
        if (~mCk2Rf(i_crk,1) && mCkJun(i_crk,1) && isempty(cHvBln_elm{i_crk,1})) || ...
           (~mCk2Rf(i_crk,2) && mCkJun(i_crk,2) && isempty(cHvBln_elm{i_crk,2}))
            
            % Hvi. bln. elements were remeshed;
            % therefore, do a more robust search
            
            tmp = [vElNew;vElHv0];
            
            [elstd,sgstd,elbln,sgbln,rmbln] = Elems_enr(mNdCrd, ...
             mLNodS(tmp,:),mCkCrd,mCkJun(i_crk,:),[0,0],he_ref);
            
            n = length(vElNew);
         
            sgstd = sgstd(elstd<=n,:);
            elstd = vElNew(elstd(elstd<=n));
            
            q = elbln{1} <= n;
            
            sgbln{1} = sgbln{1}(q,:);
            rmbln{1} = rmbln{1}(q,:);
            elbln{1} = vElNew(elbln{1}(q));
            
            q = elbln{2} <= n;
            
            sgbln{2} = sgbln{2}(q,:);
            rmbln{2} = rmbln{2}(q,:);
            elbln{2} = vElNew(elbln{2}(q));
            
        else
            
            % Hvi. enr. can only affect the new elements;
            % therefore, limit element search to 'vElNew'
            
            [elstd,sgstd,elbln,sgbln,rmbln] = Elems_enr(mNdCrd, ...
             mLNodS(vElNew,:),mCkCrd,mCkJun(i_crk,:),[0,0],he_ref);
            
            elstd    = vElNew(elstd);
            elbln{1} = vElNew(elbln{1});
            elbln{2} = vElNew(elbln{2});
            
        end
        
        % Heaviside elements to be updated
        vElHvi = [elbln{1};elstd;elbln{2}];
        mSgHvi = [sgbln{1};sgstd;sgbln{2}];
        
        % n.b. if the element to be updated is not a new element and is
        % already enriched, it is neccessary to subtract its contribution
        % from the equations if it exists
        
        % enr. el. whose eq. will need to be re-evaluated
        elupd = vElHvi(vElEnr(vElHvi) & ~vElUp0(vElHvi));
        
        if ~isempty(elupd) % subtract element enrichment from Kg and Fg
            
            % subtract Ke (se,es,ee)
            mGlStf = mGlStf - MtxStiffEnr(nGlDof,mNdCrd,mLNodS(elupd,:), ...
                1:length(elupd),cLNodE(elupd),mNDofE,vElPhz(elupd),nPhase,cDMatx, ...
                cGsEnr_omgDvS(elupd),cGsEnr_omgDvE(elupd),cGsEnr_omgWgt(elupd));
            
            % % subtract Fe (ck. face)
            % if wCkLod
            %   'Fg is to be recomputed from scratch'
            % end
            
            % subtract Fe (pre-load)
            if wPrLod
                vGlFrc = vGlFrc - ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(elupd,:),...
                    1:length(elupd),cLNodE(elupd),mNDofE,cGsEnr_omgDvS(elupd),...
                    cGsEnr_omgDvE(elupd),cGsEnr_omgWgt(elupd),mPrLod(:,elupd));
            end
            
            % subtract Fe (body)
            if wBdLod
                vGlFrc = vGlFrc - ForceEnr_Body(nGlDof,mNdCrd,mLNodS(elupd,:),...
                    1:length(elupd),cLNodE(elupd),mNDofE,cGsEnr_omgDvS(elupd),...
                    cGsEnr_omgShE(elupd),cGsEnr_omgWgt(elupd),vBdLod);
            end
            
            % clear shapes and weights
            for i = elupd(:)'
                cGsEnr_omgShS{i} = [];
                cGsEnr_omgDvS{i} = [];
                cGsEnr_omgDvE{i} = [];
                cGsEnr_omgShE{i} = [];
                cGsEnr_omgWgt{i} = [];
            end
            
            % does not work this way:
            % cGsEnr_omgShS(elupd) = [];
            
        end
        
        vElEnr(vElHvi) = true; % el. is enr.
        vElUp0(vElHvi) = true; % el. to upd.
        
        %------------------------------------------------------------------
        % Update enr. el. topology (glue new to old)
        %------------------------------------------------------------------
        
        [nNdEn8,mLNodE] = ... % get enr. topology (disjoint)
            Topo_omgEnr(nNdEn8,nEnFun_hvi,mLNodS(vElHvi,:));
        
        [nNdEn8,mLNodE] = ... % "glue" new topo. to existing
            Topo_omgEnr_add(nNdEn8,mLNodE,vElHvi,vElHv0,i_crk,0);
        
        for i = 1:length(vElHvi) % store enr. nodes (push back)
            cLNodE{vElHvi(i)}(end+1:end+nLNodE_hvi) = mLNodE(i,:);
        end
        
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        
        cNdEnr_hvi{i_crk} = ...
            [cNdEnr_hvi{i_crk};(nNdEn0+1:nNdEn8)']; nNdEn0 = nNdEn8;
        
        %------------------------------------------------------------------
        % Get bln. nodes
        %------------------------------------------------------------------
        
        q = mLNodE(ismember(vElHvi,elbln{1}),:); q = q(~rmbln{1});
        cNdBln_hvi{i_crk,1} = unique([cNdBln_hvi{i_crk,1};q(:)]);
        
        q = mLNodE(ismember(vElHvi,elbln{2}),:); q = q(~rmbln{2});
        cNdBln_hvi{i_crk,2} = unique([cNdBln_hvi{i_crk,2};q(:)]);
        
        %------------------------------------------------------------------
        % Upd. enr. data only for unique el. 
        %------------------------------------------------------------------
        
        for i = vElHvi' % (combine std. & bln.)
            cLEnDt{i}(:,end+1) = [ i_crk ; 0 ; nLNodE_hvi ; 1 ];
        end
        
        %------------------------------------------------------------------
        % Store Hvi. elements
        %------------------------------------------------------------------
        
        cHvStd_elm{i_crk}   = [cHvStd_elm{i_crk};elstd];
        cHvStd_sgm{i_crk}   = [cHvStd_sgm{i_crk};sgstd];
        
        cHvBln_elm{i_crk,1} = [elbln{1};cHvBln_elm{i_crk,1}];
        cHvBln_sgm{i_crk,1} = [sgbln{1};cHvBln_sgm{i_crk,1}];
        cHvBln_rmp{i_crk,1} = [rmbln{1};cHvBln_rmp{i_crk,1}];
        
        cHvBln_elm{i_crk,2} = [cHvBln_elm{i_crk,2};elbln{2}];
        cHvBln_sgm{i_crk,2} = [cHvBln_sgm{i_crk,2};sgbln{2}];
        cHvBln_rmp{i_crk,2} = [cHvBln_rmp{i_crk,2};rmbln{2}];
        
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        
        for i=1:length(vElHvi); ii=vElHvi(i);
            
            mElCrd = mNdCrd(mLNodS(ii,:),:);
            for j = mSgHvi(i,1):mSgHvi(i,2)
                
                [n,x] = GeoXElem(mCkCrd([j,j+1],:),mElCrd);
                
                if n == 2
                    cXsElm{ii,2}(end+1,:) = [1,2]+size(cXsElm{ii,1},1);
                end
                
                cXsElm{ii,1}(end+1:end+n,:) = x;
                
            end
        end
        
        %------------------------------------------------------------------
        
    end
end


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Safety checks
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if with_Debug

    elbad_ch1 = [];
    elbad_ch2 = [];
    elbad_ch3 = [];

    for i = find(vElEnr(:))'
        if isempty(cLNodE{i}) || isempty(cLEnDt{i})
            elbad_ch1 = [elbad_ch1;i];
        end
        if length(cLNodE{i}) ~= sum(cLEnDt{i}(3,:))
            elbad_ch2 = [elbad_ch2;i];
        end
        if size(unique(cLEnDt{i}(1:2,:)','rows'),1) ~= size(cLEnDt{i},2)
            elbad_ch3 = [elbad_ch3;i];
        end
    end

    fail = 0;

    if ~isempty(elbad_ch1)
        warning('ch1')
        fail = 1;
    end
    if ~isempty(elbad_ch2)
        warning('ch2')
        fail = 1;
    end
    if ~isempty(elbad_ch3)
        warning('ch3')
        fail = 1;
    end

    if fail
        error('safety check(s) failed')
    end

end
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


%--------------------------------------------------------------------------
% Delete "local" variables
%--------------------------------------------------------------------------

clear eltip_std elspl_std elful_std elspl_bln elful_bln elspl_rmp ...
    elful_rmp elstd sgstd elbln sgbln rmbln elupd elnew enunq 

%--------------------------------------------------------------------------
