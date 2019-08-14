
%==========================================================================
% Determine enriched shape functions for enriched elements
%==========================================================================

%--------------------------------------------------------------------------
% Initialize (already initialized)
%--------------------------------------------------------------------------
% cGsEnr_omgShE = cell(nElems,1);
% cGsEnr_omgDvE = cell(nElems,1);
%--------------------------------------------------------------------------

% upd. ck. tip. dir. (useful here)
for i = find(mCk2Up(:,1))'
    x = cCkCrd{i}(1,:)-cCkCrd{i}(2,:);
    mTpAlf(i,1) = cart2pol(x(1),x(2)); end
for i = find(mCk2Up(:,2))'
    x = cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:);
    mTpAlf(i,2) = cart2pol(x(1),x(2)); end

for i = vElUpd(:)'
    
    mElCrd = mNdCrd(mLNodS(i,:),:);
    
    % el. enr. data
    mLEnDt = cLEnDt{i};
    
    % Gauss shapes (std.)
    mGsShS = cGsEnr_omgShS{i};
    mGsDvS = cGsEnr_omgDvS{i};
    
    % gp. coord.
    mGsCrd = mGsShS*mElCrd;
    
for j = 1:size(mLEnDt,2)
    
    iCrack = mLEnDt(1,j);
    iEnTyp = mLEnDt(2,j);
    
    switch iEnTyp
            
        case 0 % Heaviside enrichment
            
            mCkCrd = cCkCrd{iCrack};
            
            if      any( i == cHvStd_elm{iCrack}    );
            
                k = i == cHvStd_elm{iCrack};
                vHvSgm = cHvStd_sgm{iCrack}(k,:);
                mHvCrd = mCkCrd(vHvSgm(1):vHvSgm(2)+1,:);
                vBlRmp = [];
                
            elseif  any( i == cHvBln_elm{iCrack,1}  );
            
                k = i == cHvBln_elm{iCrack,1};
                vHvSgm = cHvBln_sgm{iCrack,1}(k,:);
                mHvCrd = mCkCrd(vHvSgm(1):vHvSgm(2)+1,:);
                vBlRmp = [];
                
            elseif  any( i == cHvBln_elm{iCrack,2}  );
            
                k = i == cHvBln_elm{iCrack,2};
                vHvSgm = cHvBln_sgm{iCrack,2}(k,:);
                mHvCrd = mCkCrd(vHvSgm(1):vHvSgm(2)+1,:);
                vBlRmp = [];
                
            else
                error('cannot find element of enrichement type: %d, for crack: %d',[iEnTyp,iCrack])
            end
            
            % get enr. fun. values
            [mGFnVl,mNdShf] = EnrFun_Hvi(mElCrd,mGsCrd,mHvCrd);
            
            % H,x = H,y = 0
            % mGFnDv = [];
            
            if isempty(vBlRmp)
                
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtShf(mElCrd,...
                    mGsShS,mGsDvS,1,mGFnVl,[],mNdShf); % n. enr. fun = 1
                
            else
                
                vGBlVl = mGsShS*vBlRmp(:); 
                vGBlDv = mGsDvS*vBlRmp(:);
                
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtBln(mElCrd,mGsShS,...
                    mGsDvS,1,mGFnVl,[],mNdShf,vGBlVl,vGBlDv); % n. enr. fun = 1
                
            end
            
        case 1 % Branch enrichment of tip #1
            
            vTpCrd = cCkCrd{iCrack}(1,:);
            uTpAlf = mTpAlf(iCrack,1);
            
            if      any( i == cBrBln_eFl{iCrack,1} );
                
                k = i == cBrBln_eFl{iCrack,1};
                vBlRmp = cBrBln_rFl{iCrack,1}(k,:);
                  
            elseif  any( i == cBrBln_eSp{iCrack,1} );
                
                k = i == cBrBln_eSp{iCrack,1};
                vBlRmp = cBrBln_rSp{iCrack,1}(k,:);
                
            elseif  any( i == cBrStd_eFl{iCrack,1} ); % the rest can be condensed into "else"
                
                vBlRmp = [];
                
            elseif  any( i == cBrStd_eSp{iCrack,1} );
                
                vBlRmp = [];
                
            elseif  any( i == cBrStd_eTp{iCrack,1} );
                
                vBlRmp = [];
                
            else
                error('cannot find element of enrichement type: %d, for crack: %d',[iEnTyp,iCrack])
            end
            
            % get enr. fun. values
            [mGFnVl,mGFnDv,mNdShf] = EnrFun_Brn(...
                mElCrd,mGsCrd,nEnFun_brn,vTpCrd,uTpAlf);
            
            if isempty(vBlRmp)
                
                % (op1) Use jac. shape fun. to construct enr. shape fun.
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtShf(mElCrd,...
                    mGsShS,mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf);
                
            else
                
                vGBlVl = mGsShS*vBlRmp(:);
                vGBlDv = mGsDvS*vBlRmp(:);
                
                % (op1) Use jac. shape fun. to construct enr. shape fun.
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtBln(mElCrd,mGsShS,...
                    mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf,vGBlVl,vGBlDv);
                
            end
            
        case 2 % Branch enrichment of tip #2
            
            vTpCrd = cCkCrd{iCrack}(end,:);
            uTpAlf = mTpAlf(iCrack,2);
            
            if      any( i == cBrBln_eFl{iCrack,2} );
                
                k = i == cBrBln_eFl{iCrack,2};
                vBlRmp = cBrBln_rFl{iCrack,2}(k,:);
                  
            elseif  any( i == cBrBln_eSp{iCrack,2} );
                
                k = i == cBrBln_eSp{iCrack,2};
                vBlRmp = cBrBln_rSp{iCrack,2}(k,:);
                
            elseif  any( i == cBrStd_eFl{iCrack,2} ); % the rest can be condensed into "else"
                
                vBlRmp = [];
                
            elseif  any( i == cBrStd_eSp{iCrack,2} );
                
                vBlRmp = [];
                
            elseif  any( i == cBrStd_eTp{iCrack,2} );
                
                vBlRmp = [];
                
            else
                error('cannot find element of enrichement type: %d, for crack: %d',[iEnTyp,iCrack])
            end
            
            % get enr. fun. values
            [mGFnVl,mGFnDv,mNdShf] = EnrFun_Brn(...
                mElCrd,mGsCrd,nEnFun_brn,vTpCrd,uTpAlf);
            
            if isempty(vBlRmp)
                
                % (op1) Use jac. shape fun. to construct enr. shape fun.
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtShf(mElCrd,...
                    mGsShS,mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf);
                
            else
                
                vGBlVl = mGsShS*vBlRmp(:);
                vGBlDv = mGsDvS*vBlRmp(:);
                
                % (op1) Use jac. shape fun. to construct enr. shape fun.
                [nLNodE,mGsShE,mGsDvE] = ShapesEnr_WtBln(mElCrd,mGsShS,...
                    mGsDvS,nEnFun_brn,mGFnVl,mGFnDv,mNdShf,vGBlVl,vGBlDv);
                
            end
        
        otherwise
            error('unknown enrichement type: %d',iEnTyp)
    end
    
    % upd. enr. shapes
    cGsEnr_omgShE{i}(:,end+1:end+nLNodE) = mGsShE;
    cGsEnr_omgDvE{i}(:,end+1:end+nLNodE) = mGsDvE;
    
end
end