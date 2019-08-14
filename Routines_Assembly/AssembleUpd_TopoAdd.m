
%==========================================================================
% Update enriched element topology
%==========================================================================

msg1 = 'did not find split elements; i_crk = %d, i_tip = 1';
msg2 = 'did not find split elements; i_crk = %d, i_tip = 2';
msg3 = 'did not find tip elements; i_crk = %d, i_tip = 1'; % (normal)
msg4 = 'did not find tip elements; i_crk = %d, i_tip = 2'; % (normal)

% enr. el. to upd.
vElUpd = false(nElems,1);
% some were upd. already
vElUpd(vElUp0) = true;

%--------------------------------------------------------------------------
% Tip #1
%--------------------------------------------------------------------------

for i = find(mCk2Up(:,1))'
    
    % n sgm. to update
    n = mNoInc_upd(i,1);
    
    % sgm. coord. (+1 bln. sgm.)
    mCkCrd = cCkCrd{i}(1:n+2,:);
    
    if n > 0 % any enrichment
        vCkSgm = [mCkJun(i,1),n+1];
    elseif mTpRdi(i,1) > 0 % Brn. enr.
        vCkSgm = [0,0];
    else % no enr.
        mCkCrd = [];
        vCkSgm = [];
    end
    
    [vHvStd_elm,mHvStd_sgm,mHvBln_elm,mHvBln_sgm,mHvBln_rmp,mBrStd_eTp,...
     mBrStd_eSp,mBrStd_eFl,mBrBln_eSp,mBrBln_eFl,mBrBln_rSp,mBrBln_rFl]=...
     Elems_enr(mNdCrd,mLNodS,mCkCrd,vCkSgm,[mTpRdi(i,1),0],he_ref);
    
    if ~isempty(mBrStd_eTp{1})
        
        vElBrn = [mBrStd_eTp{1}; mBrStd_eSp{1}; ...
            mBrStd_eFl{1}; mBrBln_eSp{1}; mBrBln_eFl{1}];
        
        % el. whose Ke to negate
        p = vElBrn(vElEnr(vElBrn) & ~vElUpd(vElBrn));
        
        if ~isempty(p)
        
            %--------------------------------------------------------------
            % Negate Ke from Kg
            %--------------------------------------------------------------
            
            mGlStf = mGlStf - ...
                MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));
            
            %--------------------------------------------------------------
            % Negate Fe from Fg
            %--------------------------------------------------------------
            if wPrLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));
                
            end
            if wBdLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);
                
            end
            %--------------------------------------------------------------
            % Negate shapes
            %--------------------------------------------------------------
            for k = p(:)'
                
                % remove shapes (std.)
                cGsEnr_omgShS{k} = [];
                cGsEnr_omgDvS{k} = [];
                % remove shapes (enr.)
                cGsEnr_omgDvE{k} = [];
                cGsEnr_omgShE{k} = [];
                % remove weights
                cGsEnr_omgWgt{k} = [];
                
            end
            %--------------------------------------------------------------
            
        end
        
        vElEnr(vElBrn) = true;
        vElUpd(vElBrn) = true;
        
        %------------------------------------------------------------------
        % Update enr. el. topology
        %------------------------------------------------------------------
        
        [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,nEnFun_brn,mLNodS(vElBrn,:));
        
        for j = 1:length(vElBrn)
            cLNodE{vElBrn(j)}(end+1:end+nLNodE_brn) = mLNodE(j,:);
        end
        
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        
        cNdEnr_brn{i,1} = (nNdEn0+1:nNdEn8)'; nNdEn0 = nNdEn8;
        
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
            
            [n,x] = GeoXElem(mCkCrd([1,2],:), ...
                mNdCrd(mLNodS(j,:),:));
            
            if n == 2
                cXsElm{j,2}(end+1,:) = [1,2]+size(cXsElm{j,1},1);
            end
            
            cXsElm{j,1}(end+1:end+n,:) = x;
            
        end
        
        %------------------------------------------------------------------
        
    else
        mTpRdi(i,1) = 0;
        mTpAct(i,1) = 0;
        % warning(msg3,i)
    end
    
    if mNoInc_upd(i,1) > 0 && ~isempty(vHvStd_elm)
        
        % old el.
        vElHv0 = [cHvBln_elm{i,1};cHvStd_elm{i}];
        mSgHv0 = [cHvBln_sgm{i,1};cHvStd_sgm{i}];
        
        % new el. (overlap with vElHv0 - to extra safe)
        vElHvi = [mHvBln_elm{1};vHvStd_elm;mHvBln_elm{2}];
        mSgHvi = [mHvBln_sgm{1};mHvStd_sgm;mHvBln_sgm{2}];
        elunq  = vElHvi(~ismember(vElHvi,vElHv0)); % unq.
        
        if isempty(elunq)
            warning('Could not find new (unique) Hvi. elements');
        end
        
        % el. whose Ke to negate (safe but overkill)
        p = vElHvi(vElEnr(vElHvi) & ~vElUpd(vElHvi));
        
        if ~isempty(p)
            
            %--------------------------------------------------------------
            % Negate Ke from Kg
            %--------------------------------------------------------------
            
            mGlStf = mGlStf - ...
                MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));
            
            %--------------------------------------------------------------
            % Negate Fe from Fg
            %--------------------------------------------------------------
            if wPrLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));
                
            end
            if wBdLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);
                
            end
            %--------------------------------------------------------------
            % Negate shapes
            %--------------------------------------------------------------
            for k = p(:)'
                
                % remove shapes (std.)
                cGsEnr_omgShS{k} = [];
                cGsEnr_omgDvS{k} = [];
                % remove shapes (enr.)
                cGsEnr_omgDvE{k} = [];
                cGsEnr_omgShE{k} = [];
                % remove weights
                cGsEnr_omgWgt{k} = [];
                
            end
            %--------------------------------------------------------------
            
        end
        
        vElEnr([elunq;p]) = true;
        vElUpd([elunq;p]) = true;
        
        if ~isempty(elunq)
            
            %--------------------------------------------------------------
            % Update enr. el. topology
            %--------------------------------------------------------------
            %{
            % old topo. (std. for ref.)
            tmp1 = mLNodS(vElHv0,:);
            % new topo. (std. for ref.)
            tmp2 = mLNodS(vElHvi,:);

            % get nodes that are shared
            [p_mbr,p_ref] = ismember(tmp2(:),tmp1(:));

            % ind. into tmp1
            p_ref = p_ref(p_mbr);
            % ind. into tmp2
            p_mbr = find(p_mbr);

            % ind -> sub into tmp1
            iLNodS = ceil(p_ref/nElHv0);
            iElHv0 = p_ref-(iLNodS-1)*nElHv0;

            % get new enriched element topology
            [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,1,tmp2);

            % bln. nd. must be in new enr. topo.
            vNdBln = unique(mLNodE(p_mbr));
            nNdBln = length(vNdBln);

            % rmv. excess enr. nd.
            nNdEn8 = nNdEn8-nNdBln;

            % shift nd. down since nd. rmv.
            for j = nNdBln:-1:1
                p = mLNodE(:) > vNdBln(j);
                mLNodE(p) = mLNodE(p) - 1;
            end

            % make topo. consistent
            for j = 1:length(p_mbr)

                % old element
                k = vElHv0(iElHv0(j));
                % old element enr. data
                mLEnDt = cLEnDt{k};

                % enr. nd. position inside enr. el. nodes
                p = find(mLEnDt(1,:)==i & mLEnDt(2,:)==0);
                p = sum(mLEnDt(3,1:p-1)) + iLNodS(j);

                % modify to be consistent
                mLNodE(p_mbr(j)) = cLNodE{k}(p);

            end

            for j = jElUnq(:)' % only add if el. is not in vElHv0
                cLNodE{vElHvi(j)}(end+1:end+nLNodS) = mLNodE(j,:);
            end
            %}
            
            [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,nEnFun_hvi,mLNodS(elunq,:));
            [nNdEn8,mLNodE] = Topo_omgEnr_add(nNdEn8,mLNodE,elunq,vElHv0,i,0);
            
            for j = 1:length(elunq) % nLNodS == nLNodE_hvi
                cLNodE{elunq(j)}(end+1:end+nLNodS) = mLNodE(j,:);
            end
            
            %--------------------------------------------------------------
            % Get enr. nodes
            %--------------------------------------------------------------
            
            cNdEnr_hvi{i} = [cNdEnr_hvi{i};(nNdEn0+1:nNdEn8)']; nNdEn0=nNdEn8;
            
            %--------------------------------------------------------------
            % Get bln. nodes
            %--------------------------------------------------------------
            
            tmp = mLNodE(ismember(elunq,mHvBln_elm{1}),:);
            tmp = unique(tmp(~mHvBln_rmp{1}));
            cNdBln_hvi{i,1} = tmp(:);
            
            %--------------------------------------------------------------
            % Upd. enr. data only for unique el. (combine std. & bln.)
            %--------------------------------------------------------------
            
            for j = elunq(:)' % only add if el. is not in vElHv0
                cLEnDt{j}(:,end+1) = [ i ; 0 ; nLNodS ; 1 ];
            end
            
            %--------------------------------------------------------------
        
        end
        
        %------------------------------------------------------------------
        % Store Hvi. enrichment data
        %------------------------------------------------------------------
        
        p = ~ismember(vElHv0,vElHvi);
        
        vElHv0 = vElHv0(p);
        mSgHv0 = mSgHv0(p,:);
        
        % replace with new
        cHvBln_elm{i,1} = mHvBln_elm{1};
        cHvBln_sgm{i,1} = mHvBln_sgm{1};
        cHvBln_rmp{i,1} = mHvBln_rmp{1};
        
        % push forward
        cHvStd_elm{i} = [vHvStd_elm; mHvBln_elm{2} ; vElHv0];
        cHvStd_sgm{i} = [mHvStd_sgm; mHvBln_sgm{2} ; mSgHv0];
        
        % upd. old
        cHvBln_sgm{i,2} = cHvBln_sgm{i,2};
        
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        
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
        
        %------------------------------------------------------------------
        
    else
        if mNoInc_upd(i,1) > 0
            warning(msg1,i)
        end
    end
end

%--------------------------------------------------------------------------
% Tip #2
%--------------------------------------------------------------------------

for i = find(mCk2Up(:,2))'
    
    % n sgm. to update
    n = mNoInc_upd(i,2);
    
    % ((size(cCkCrd{i},1)-1)-n)-1;
    i_sg0  = size(cCkCrd{i},1)-n-2;     % reference segment
    mCkCrd = cCkCrd{i}(i_sg0+1:end,:);  % +1 blending segment
    
    if n > 0 % any enrichment
        if mCkJun(i,2) == 0; vCkSgm = [1,0];
        else vCkSgm = [1,mCkJun(i,2)-i_sg0]; end
    elseif mTpRdi(i,2) > 0 % Brn. enr.
        vCkSgm = [0,0];
    else % no enr.
        mCkCrd = [];
        vCkSgm = [];
    end
    
    [vHvStd_elm,mHvStd_sgm,mHvBln_elm,mHvBln_sgm,mHvBln_rmp,mBrStd_eTp,...
     mBrStd_eSp,mBrStd_eFl,mBrBln_eSp,mBrBln_eFl,mBrBln_rSp,mBrBln_rFl]=...
     Elems_enr(mNdCrd,mLNodS,mCkCrd,vCkSgm,[0,mTpRdi(i,2)],he_ref);
    
    if ~isempty(mBrStd_eTp{2})
        
        vElBrn = [mBrStd_eTp{2}; mBrStd_eSp{2}; ...
            mBrStd_eFl{2}; mBrBln_eSp{2}; mBrBln_eFl{2}];
        
        % el. whose Ke to negate
        p = vElBrn(vElEnr(vElBrn) & ~vElUpd(vElBrn));
        
        if ~isempty(p)
            
            %--------------------------------------------------------------
            % Negate Ke from Kg
            %--------------------------------------------------------------
            
            mGlStf = mGlStf - ...
                MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));
            
            %--------------------------------------------------------------
            % Negate Fe from Fg
            %--------------------------------------------------------------
            if wPrLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));
                
            end
            if wBdLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);
                
            end
            %--------------------------------------------------------------
            % Negate shapes
            %--------------------------------------------------------------
            for k = p(:)'
                
                % remove shapes (std.)
                cGsEnr_omgShS{k} = [];
                cGsEnr_omgDvS{k} = [];
                % remove shapes (enr.)
                cGsEnr_omgDvE{k} = [];
                cGsEnr_omgShE{k} = [];
                % remove weights
                cGsEnr_omgWgt{k} = [];
                
            end
            %--------------------------------------------------------------
            
        end
        
        vElEnr(vElBrn) = true;
        vElUpd(vElBrn) = true;
        
        %------------------------------------------------------------------
        % Update enr. el. topology
        %------------------------------------------------------------------
        
        [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,nEnFun_brn,mLNodS(vElBrn,:));
        
        for j = 1:length(vElBrn)
            cLNodE{vElBrn(j)}(end+1:end+nLNodE_brn) = mLNodE(j,:);
        end
        
        %------------------------------------------------------------------
        % Get enr. nodes
        %------------------------------------------------------------------
        
        cNdEnr_brn{i,2} = (nNdEn0+1:nNdEn8)'; nNdEn0 = nNdEn8;
        
        %------------------------------------------------------------------
        % Save el. enr. data: ck. id, enr. id, n. enr. node, GP.
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
            
            [n,x] = GeoXElem(mCkCrd([end,end-1],:), ...
                mNdCrd(mLNodS(j,:),:));
            
            if n == 2
                cXsElm{j,2}(end+1,:) = [1,2]+size(cXsElm{j,1},1);
            end
            
            cXsElm{j,1}(end+1:end+n,:) = x;
            
        end
        
        %------------------------------------------------------------------
        
    else
        mTpRdi(i,2) = 0;
        mTpAct(i,2) = 0;
        % warning(msg4,i)
    end
    
    if mNoInc_upd(i,2) > 0 && ~isempty(vHvStd_elm)
        
        % old el.
        vElHv0 = [cHvStd_elm{i};cHvBln_elm{i,2}];
        mSgHv0 = [cHvStd_sgm{i};cHvBln_sgm{i,2}];
        
        % new el. (overlap with vElHv0 - to extra safe)
        vElHvi = [mHvBln_elm{1};vHvStd_elm;mHvBln_elm{2}];
        mSgHvi = [mHvBln_sgm{1};mHvStd_sgm;mHvBln_sgm{2}];
        elunq = vElHvi(~ismember(vElHvi,vElHv0)); % unq.
        
        if isempty(elunq)
            warning('Could not find new (unique) Hvi. elements');
        end
        
        % el. whose Ke to negate (safe but overkill)
        p = vElHvi(vElEnr(vElHvi) & ~vElUpd(vElHvi));
        
        % negate Ke from Kg.
        if ~isempty(p)
            
            %--------------------------------------------------------------
            % Negate Ke from Kg
            %--------------------------------------------------------------
            
            mGlStf = mGlStf - ...
                MtxStiffEnr(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                vElPhz(p),nPhase,cDMatx,cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p));
            
            %--------------------------------------------------------------
            % Negate Fe from Fg
            %--------------------------------------------------------------
            if wPrLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Resid(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgDvE(p),cGsEnr_omgWgt(p),mPrLod(:,p));
                
            end
            if wBdLod
                
                vGlFrc = vGlFrc - ...
                    ForceEnr_Body(nGlDof,mNdCrd,mLNodS(p,:),1:length(p),cLNodE(p),mNDofE,...
                    cGsEnr_omgDvS(p),cGsEnr_omgShE(p),cGsEnr_omgWgt(p),vBdLod);
                
            end
            %--------------------------------------------------------------
            % Negate shapes
            %--------------------------------------------------------------
            for k = p(:)'
                
                % remove shapes (std.)
                cGsEnr_omgShS{k} = [];
                cGsEnr_omgDvS{k} = [];
                % remove shapes (enr.)
                cGsEnr_omgDvE{k} = [];
                cGsEnr_omgShE{k} = [];
                % remove weights
                cGsEnr_omgWgt{k} = [];
                
            end
            %--------------------------------------------------------------
            
        end
        
        vElEnr([elunq;p]) = true;
        vElUpd([elunq;p]) = true;
        
        if ~isempty(elunq)
            
            %--------------------------------------------------------------
            % Update enr. el. topology
            %--------------------------------------------------------------
            %{
            % old topo. (std. for ref.)
            tmp1 = mLNodS(vElHv0,:);
            % new topo. (std. for ref.)
            tmp2 = mLNodS(vElHvi,:);

            % get nodes that are shared
            [p_mbr,p_ref] = ismember(tmp2(:),tmp1(:));

            % ind. into tmp1
            p_ref = p_ref(p_mbr);
            % ind. into tmp2
            p_mbr = find(p_mbr);

            % ind -> sub into tmp1
            iLNodS = ceil(p_ref/nElHv0);
            iElHv0 = p_ref-(iLNodS-1)*nElHv0;

            % get new enriched element topology
            [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,1,tmp2);

            % bln. nd. must be in new enr. topo.
            vNdBln = unique(mLNodE(p_mbr));
            nNdBln = length(vNdBln);

            % rmv. excess enr. nd.
            nNdEn8 = nNdEn8-nNdBln;

            % shift nd. down since nd. rmv.
            for j = nNdBln:-1:1
                p = mLNodE(:) > vNdBln(j);
                mLNodE(p) = mLNodE(p) - 1;
            end

            % make topo. consistent
            for j = 1:length(p_mbr)

                % old element
                k = vElHv0(iElHv0(j));
                % old element enr. data
                mLEnDt = cLEnDt{k};

                % enr. position
                p = find(mLEnDt(1,:)==i & mLEnDt(2,:)==0);
                p = sum(mLEnDt(3,1:p-1)) + iLNodS(j);

                % modify to be consistent
                mLNodE(p_mbr(j)) = cLNodE{k}(p);

            end

            for j = jElUnq(:)' % only add if el. is not in vElHv0
                cLNodE{vElHvi(j)}(end+1:end+nLNodS) = mLNodE(j,:);
            end
            %}
            
            [nNdEn8,mLNodE] = Topo_omgEnr(nNdEn8,nEnFun_hvi,mLNodS(elunq,:));
            [nNdEn8,mLNodE] = Topo_omgEnr_add(nNdEn8,mLNodE,elunq,vElHv0,i,0);
            
            for j = 1:length(elunq) % nLNodS == nLNodE_hvi
                cLNodE{elunq(j)}(end+1:end+nLNodS) = mLNodE(j,:);
            end
            
            %--------------------------------------------------------------
            % Get enr. nodes
            %--------------------------------------------------------------
            
            cNdEnr_hvi{i} = [cNdEnr_hvi{i};(nNdEn0+1:nNdEn8)']; nNdEn0=nNdEn8;
            
            %--------------------------------------------------------------
            % Get bln. nodes
            %--------------------------------------------------------------
            
            tmp = mLNodE(ismember(elunq,mHvBln_elm{2}),:);
            tmp = unique(tmp(~mHvBln_rmp{2}));
            cNdBln_hvi{i,2} = tmp(:);
            
            %--------------------------------------------------------------
            % Upd. enr. data only for unique el. (combine std. & bln.)
            %--------------------------------------------------------------
            
            for j = elunq(:)' % only add if el. is not in vElHv0
                cLEnDt{j}(:,end+1) = [ i ; 0 ; nLNodS ; 1 ];
            end
            
            %--------------------------------------------------------------
            
        end
        
        %------------------------------------------------------------------
        % Store Hvi. enrichment data
        %------------------------------------------------------------------
        
        p = ~ismember(vElHv0,vElHvi);
        
        vElHv0 = vElHv0(p);
        mSgHv0 = mSgHv0(p,:);
        
        % replace with new
        cHvBln_elm{i,2} = mHvBln_elm{2};
        cHvBln_sgm{i,2} = mHvBln_sgm{2}+i_sg0;
        cHvBln_rmp{i,2} = mHvBln_rmp{2};
        
        % push back
        cHvStd_elm{i} = [vElHv0; mHvBln_elm{1}        ; vHvStd_elm     ];
        cHvStd_sgm{i} = [mSgHv0; mHvBln_sgm{1}+i_sg0 ; mHvStd_sgm+i_sg0];
        
        %------------------------------------------------------------------
        % Determine background mesh
        %------------------------------------------------------------------
        
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
        
        %------------------------------------------------------------------
        
    else
        if mNoInc_upd(i,2) > 0
            warning(msg2,i)
        end
    end
    
end


% get std. el.
vElStd = find(~vElEnr)';
nElStd = length(vElStd);

% get enr. el.
vElEnr = find(vElEnr)';
nElEnr = length(vElEnr);

% get enr. el. to upd.
vElUpd = find(vElUpd);
nElUpd = length(vElUpd);

% cleanup sub-cells
for i = vElUpd(:)'
    if ~isempty(cXsElm{i,1})
        
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
end
