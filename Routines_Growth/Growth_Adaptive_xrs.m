
%--------------------------------------------------------------------------
% Adaptive fracture advance/remeshing based on intersection proximity
%--------------------------------------------------------------------------

% mTpRdi == 0, implies proximate intersection
% remeshing is allowed even if mNoInc_grw == 0

mCk2Rf = mCk2Up & mTpRdi == 0 & mCkRfN<nRefine_xrs; % (mCk2Rf is local)
if any(mCk2Rf(:)); fprintf('\n\tadaptive remeshing: pre-intersection\n')
    
    % assert tips to be refined
    mCkRfn(mCk2Rf) = true;
    % update n. times to refine
    mCkRfN = mCkRfN+real(mCk2Rf);
    
    % set new enrichment radii: halfing
    mTpRdi(mCk2Rf) = mTpRdi_ref(mCk2Rf);
    mTpRdi = mTpRdi.*(0.5.^real(mCk2Rf));
    
    % restore previous tip activities
    mTpAct(mCk2Rf) = mTpAct_ref(mCk2Rf);
    
    % nullify junctions
    mCkJun(mCk2Rf) = 0;
    
    % go to crack before intersection
    for i = find(mCk2Rf(:,2) & mNoInc_grw(:,2)>0)'
        
        % pulling back and halfing the crack increment
        cCkCrd{i} = cCkCrd{i}(1:end-mNoInc_grw(i,2)+1,:);
        d = cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:); d=d./norm(d);
        cCkCrd{i}(end,:) = cCkCrd{i}(end-1,:)+(f_inc*mTpRdi(i,2))*d;
        
        % growth by 1 inc.
        mNoInc_grw(i,2) = 1;
        
    end
    
    % go to crack before intersection
    for i = find(mCk2Rf(:,1) & mNoInc_grw(:,1)>0)'
        
        % pulling back and halfing the crack increment
        cCkCrd{i} = cCkCrd{i}(mNoInc_grw(i,1):end,:);
        d = cCkCrd{i}(1,:)-cCkCrd{i}(2,:); d=d./norm(d);
        cCkCrd{i}(1,:) = cCkCrd{i}(2,:)+(f_inc*mTpRdi(i,1))*d;
        
        if mCkJun(i,2) ~= 0
            mCkJun(i,2) = mCkJun(i,2)-mNoInc_grw(i,1)+1;
        end
        
        % growth by 1 inc.
        mNoInc_grw(i,1) = 1;
        
    end
    
    % at least 1 inc. to upd.
    mNoInc_upd(mCk2Rf) = 1;
    
end