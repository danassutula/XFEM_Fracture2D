
%==========================================================================
% Determine enriched element topology
%   
%   enrichment sequence:
%       (1) - Branch, tip #1
%       (2) - Branch, tip #2
%       (3) - Heaviside
%
%==========================================================================

%--------------------------------------------------------------------------
% Initialize
%--------------------------------------------------------------------------
% enr. el.
vElEnr = false(1,nElems);

% is ck. discretized
vCkRmv = false(nCrack,1);

% Hvi. enrichment
cHvStd_elm = cell(nCrack,1);
cHvStd_sgm = cell(nCrack,1);
cHvBln_elm = cell(nCrack,2);
cHvBln_sgm = cell(nCrack,2);
cHvBln_rmp = cell(nCrack,2);

% Brn. enrichment
cBrStd_eTp = cell(nCrack,2);
cBrStd_eSp = cell(nCrack,2);
cBrStd_eFl = cell(nCrack,2);
cBrBln_eSp = cell(nCrack,2);
cBrBln_eFl = cell(nCrack,2);
cBrBln_rSp = cell(nCrack,2);
cBrBln_rFl = cell(nCrack,2);

% enr. el. topology
cLNodE = cell(nElems,1); % {[q_enrNodes];...}
% enr. track
cLEnDt = cell(nElems,1); % {[i_crk,i_enr,n_enrNodes,id_quadRule]',...}
% for el. subdiv.
cXsElm = cell(nElems,2); % {[intersection points]},{[constrained edges]}

% enr. nodes
cNdEnr_hvi = cell(nCrack,1);
cNdEnr_brn = cell(nCrack,2);

% bln. nodes
cNdBln_hvi = cell(nCrack,1);
cNdBln_brn = cell(nCrack,2); % not used

% enr. node count
nNdEnr = 0;
nNdEn0 = 0;

%--------------------------------------------------------------------------
% Assemble el. enr. data
%--------------------------------------------------------------------------
for i = 1:nCrack

mCkCrd = cCkCrd{i};

[vHvStd_elm,mHvStd_sgm,mHvBln_elm,mHvBln_sgm,mHvBln_rmp,mBrStd_eTp,...
 mBrStd_eSp,mBrStd_eFl,mBrBln_eSp,mBrBln_eFl,mBrBln_rSp,mBrBln_rFl] = ...
 Elems_enr(mNdCrd,mLNodS,mCkCrd,mCkJun(i,:),mTpRdi(i,:),he_ref);

if ~isempty(vHvStd_elm)
    % need: for j = 1:2
    if ~isempty(mBrStd_eTp{1})
        
        %------------------------------------------------------------------
        % Get enr. el. topology
        %------------------------------------------------------------------
        vElBrn = [mBrStd_eTp{1}; mBrStd_eSp{1}; ...
            mBrStd_eFl{1}; mBrBln_eSp{1}; mBrBln_eFl{1}];
        
        % el. is enr.
        vElEnr(vElBrn) = true;
        
        [nNdEnr,mLNodE] = ...
            Topo_omgEnr(nNdEnr,nEnFun_brn,mLNodS(vElBrn,:));
        
        for j = 1:length(vElBrn)
            cLNodE{vElBrn(j)}(end+1:end+nLNodE_brn) = mLNodE(j,:);
        end
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        cNdEnr_brn{i,1} = (nNdEn0+1:nNdEnr)'; nNdEn0 = nNdEnr;
        %------------------------------------------------------------------
        % Get bln. nodes
        %------------------------------------------------------------------
%         p = [mBrBln_eSp{1};mBrBln_eFl{1}]; % bln. el.
%         q = [mBrBln_rSp{1};mBrBln_rFl{1}]; % bln. ramp
%         
%         % get bln. el. topo. (input is sorted)
%         mLNodE = mLNodE(ismember(vElBrn,p),:);
%         
%         % bln. node id
%         k = find(~q(:));
%         % get bln. nodes for iEnFun = 1
%         [~,j] = unique(mLNodE(k)); k = k(j);
%         
%         % n. nd. in mLNodE per enrichment
%         n = length(p)*nLNodS;
%         % n. bln. nodes
%         m = length(k);
%         
%         for j = 1:nEnFun_brn
%             cNdBln_brn{i,1}(end+1:end+m,1) = mLNodE(k); k = k + n;
%         end
        %------------------------------------------------------------------
        % Save el. enr. data: ck. id, enr. id, n. enr. node, GP.
        %------------------------------------------------------------------
        for j = mBrStd_eTp{1}'
            cLEnDt{j}(:,end+1) = [ i ; 1 ; nLNodE_brn ; 3]; end
        for j = mBrStd_eSp{1}'
            cLEnDt{j}(:,end+1) = [ i ; 1 ; nLNodE_brn ; 2]; end
        for j = mBrStd_eFl{1}'
            cLEnDt{j}(:,end+1) = [ i ; 1 ; nLNodE_brn ; 2]; end
        for j = mBrBln_eSp{1}'
            cLEnDt{j}(:,end+1) = [ i ; 1 ; nLNodE_brn ; 2]; end
        for j = mBrBln_eFl{1}'
            cLEnDt{j}(:,end+1) = [ i ; 1 ; nLNodE_brn ; 2]; end
        %------------------------------------------------------------------
        % Save Brn. enrichment data: el., bln. ramp
        %------------------------------------------------------------------
        cBrStd_eTp{i,1} = mBrStd_eTp{1};
        cBrStd_eSp{i,1} = mBrStd_eSp{1};
        cBrStd_eFl{i,1} = mBrStd_eFl{1};
        
        cBrBln_eSp{i,1} = mBrBln_eSp{1};
        cBrBln_rSp{i,1} = mBrBln_rSp{1};
        cBrBln_eFl{i,1} = mBrBln_eFl{1};
        cBrBln_rFl{i,1} = mBrBln_rFl{1};
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        for j = mBrStd_eTp{1}'
            
            [n,x] = GeoXElem(mCkCrd([1,2],:),mNdCrd(mLNodS(j,:),:));
            
            if n == 2
                cXsElm{j,2}(end+1,:) = [1,2]+size(cXsElm{j,1},1);
            end
            
            cXsElm{j,1}(end+1:end+n,:) = x;
            
        end
        %------------------------------------------------------------------
        
    else
        mTpRdi(i,1) = 0;
        mTpAct(i,1) = 0;
    end
    
    if ~isempty(mBrStd_eTp{2})
        
        %------------------------------------------------------------------
        % Get enr. el. topology
        %------------------------------------------------------------------
        vElBrn = [mBrStd_eTp{2}; mBrStd_eSp{2}; ...
            mBrStd_eFl{2}; mBrBln_eSp{2}; mBrBln_eFl{2}];
        
        % el. is enr.
        vElEnr(vElBrn) = true;
        
        [nNdEnr,mLNodE] = ...
            Topo_omgEnr(nNdEnr,nEnFun_brn,mLNodS(vElBrn,:));
        
        for j = 1:length(vElBrn)
            cLNodE{vElBrn(j)}(end+1:end+nLNodE_brn) = mLNodE(j,:);
        end
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        cNdEnr_brn{i,2} = (nNdEn0+1:nNdEnr)'; nNdEn0 = nNdEnr;
        %------------------------------------------------------------------
        % Get bln. nodes
        %------------------------------------------------------------------
%         p = [mBrBln_eSp{2};mBrBln_eFl{2}]; % bln. el.
%         q = [mBrBln_rSp{2};mBrBln_rFl{2}]; % bln. ramp
%         
%         % get bln. el. topo. (input is sorted)
%         mLNodE = mLNodE(ismember(vElBrn,p),:);
%         
%         % bln. node id
%         k = find(~q(:));
%         % get bln. nodes for iEnFun = 1
%         [~,j] = unique(mLNodE(k)); k = k(j);
%         
%         % n. nd. in mLNodE per enrichment
%         n = length(p)*nLNodS;
%         % n. bln. nodes
%         m = length(k);
%         
%         for j = 1:nEnFun_brn
%             cNdBln_brn{i,2}(end+1:end+m,1) = mLNodE(k); k = k + n;
%         end
        %------------------------------------------------------------------
        % Save el. enrichment data: ck. id, enr. id, n. enr. node, GP.
        %------------------------------------------------------------------
        for j = mBrStd_eTp{2}'
            cLEnDt{j}(:,end+1) = [ i ; 2 ; nLNodE_brn ; 3]; end
        for j = mBrStd_eSp{2}'
            cLEnDt{j}(:,end+1) = [ i ; 2 ; nLNodE_brn ; 2]; end
        for j = mBrStd_eFl{2}'
            cLEnDt{j}(:,end+1) = [ i ; 2 ; nLNodE_brn ; 2]; end
        for j = mBrBln_eSp{2}'
            cLEnDt{j}(:,end+1) = [ i ; 2 ; nLNodE_brn ; 2]; end
        for j = mBrBln_eFl{2}'
            cLEnDt{j}(:,end+1) = [ i ; 2 ; nLNodE_brn ; 2]; end
        %------------------------------------------------------------------
        % Save Brn. enrichment data: el., bln. ramp
        %------------------------------------------------------------------    
        cBrStd_eTp{i,2} = mBrStd_eTp{2};
        cBrStd_eSp{i,2} = mBrStd_eSp{2};
        cBrStd_eFl{i,2} = mBrStd_eFl{2};

        cBrBln_eSp{i,2} = mBrBln_eSp{2};
        cBrBln_rSp{i,2} = mBrBln_rSp{2};
        cBrBln_eFl{i,2} = mBrBln_eFl{2};
        cBrBln_rFl{i,2} = mBrBln_rFl{2};
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        for j = mBrStd_eTp{2}'
            
            [n,x] = GeoXElem(mCkCrd([end,end-1],:),mNdCrd(mLNodS(j,:),:));
            
            if n == 2
                cXsElm{j,2}(end+1,:) = [1,2]+size(cXsElm{j,1},1);
            end
            
            cXsElm{j,1}(end+1:end+n,:) = x;
            
        end
        %------------------------------------------------------------------
        
    else
        mTpRdi(i,2) = 0;
        mTpAct(i,2) = 0;
    end
    
    %----------------------------------------------------------------------
    % Get enr. el. topology
    %----------------------------------------------------------------------
    vElHvi = [mHvBln_elm{1};vHvStd_elm;mHvBln_elm{2}];
    mSgHvi = [mHvBln_sgm{1};mHvStd_sgm;mHvBln_sgm{2}]; % for later
    
    % el. is enr.
    vElEnr(vElHvi) = true;
    
    [nNdEnr,mLNodE] = Topo_omgEnr(nNdEnr,1,mLNodS(vElHvi,:));
    
    for j = 1:length(vElHvi)
        cLNodE{vElHvi(j)}(end+1:end+nLNodS) = mLNodE(j,:);
    end
    %----------------------------------------------------------------------
    % Get enr. nodes
    %----------------------------------------------------------------------
    cNdEnr_hvi{i} = (nNdEn0+1:nNdEnr)'; nNdEn0 = nNdEnr;
    %----------------------------------------------------------------------
    % Get bln. nodes
    %----------------------------------------------------------------------
    % end #1
    p = mLNodE(ismember(vElHvi,mHvBln_elm{1}),:);
    p = unique(p(~mHvBln_rmp{1}));
    cNdBln_hvi{i,1} = p(:);
    % end #2
    p = mLNodE(ismember(vElHvi,mHvBln_elm{2}),:);
    p = unique(p(~mHvBln_rmp{2}));
    cNdBln_hvi{i,2} = p(:);
    %----------------------------------------------------------------------
    % Get el. enrichment data: ck. id, enr. id, n. enr. node, GP.
    %----------------------------------------------------------------------
    for j = mHvBln_elm{1}'
        cLEnDt{j}(:,end+1) = [ i ; 0 ; nLNodS ; 1]; end
    for j = vHvStd_elm'
        cLEnDt{j}(:,end+1) = [ i ; 0 ; nLNodS ; 1]; end
    for j = mHvBln_elm{2}'
        cLEnDt{j}(:,end+1) = [ i ; 0 ; nLNodS ; 1]; end
    %----------------------------------------------------------------------
    % Save Brn. enrichment data: el., bln. ramp, ck. sgm.
    %----------------------------------------------------------------------
    cHvStd_elm{i}   = vHvStd_elm;
    cHvStd_sgm{i}   = mHvStd_sgm;
    
    cHvBln_elm{i,1} = mHvBln_elm{1};
    cHvBln_sgm{i,1} = mHvBln_sgm{1};
    cHvBln_rmp{i,1} = mHvBln_rmp{1};
    
    cHvBln_elm{i,2} = mHvBln_elm{2};
    cHvBln_sgm{i,2} = mHvBln_sgm{2};
    cHvBln_rmp{i,2} = mHvBln_rmp{2};
    %----------------------------------------------------------------------
    % Determine background mesh
    %----------------------------------------------------------------------
    for j=1:length(vElHvi); jj=vElHvi(j);
        
        mElCrd = mNdCrd(mLNodS(jj,:),:);
        for k = mSgHvi(j,1):mSgHvi(j,2)
            
            [n,x] = GeoXElem(mCkCrd([k,k+1],:),mElCrd);
            
            if n == 2
                cXsElm{jj,2}(end+1,:) = [1,2]+size(cXsElm{jj,1},1);
            end
            
            cXsElm{jj,1}(end+1:end+n,:) = x;
            
        end
    end
    %----------------------------------------------------------------------
    
else
    warning('mesh is too coarse; did not find i_crk = %d',i)
    mTpRdi(i,:) = 0; mTpAct(i,:) = 0; vCkRmv(i) = true;
end
end


% get std. el.
vElStd = find(~vElEnr);
nElStd = length(vElStd);

% get enr. el.
vElEnr = find(vElEnr);
nElEnr = length(vElEnr);

% undiscretized cracks
vCkRmv = find(vCkRmv);
nCkRmv = length(vCkRmv);


% cleanup sub-cells
for i = unique([cell2mat(cBrStd_eTp(:)); ...
        cell2mat(cHvBln_elm(:));cell2mat(cHvStd_elm(:))])';
    
    [q_unq,q_rpt,q_r2u] = Nodes_unq(cXsElm{i,1},tol_abs);
    % unique points, repeated points, q_unq=q_r2u(q_rpt)
    
    if ~isempty(q_rpt) % found repeating points
        
        for j = q_rpt(:)'
            % replace repeated points with unique
            cXsElm{i,2}(cXsElm{i,2}==j) = q_r2u(j);
        end
        
        for j = find(1:length(q_unq) ~= q_unq(:)')
            % arrange points consequently: 1,2,...
            cXsElm{i,2}(cXsElm{i,2}==q_unq(j)) = j;
        end
        
        % remove repeated constrained edges
        cXsElm{i,2} = unique(sort(cXsElm{i,2},2),'rows');
        cXsElm{i,1} = cXsElm{i,1}(q_unq,:);
        
    end
end
