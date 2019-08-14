
%--------------------------------------------------------------------------
% Fracture advance/remeshing based on the maximum number of refinements
%--------------------------------------------------------------------------

if nRefine_inc ~= 0 % (do only if possible to refine)

% cracks incremented for the first time
mCk2Rf = ~mIsInc & mNoInc_grw>0 & mCkRfN<nRefine_inc;

% n.b. mIsInc==0 says that a crack either has not propagated yet or that 
% is has just propagated; i.e. has grown for the first time

if any(mCk2Rf(:))
    
    % update tip positions: halfing increment max times
    fprintf('\n\tadaptive remeshing: refining (max)\n')
    
    % first pull back intersections
    q = mCk2Rf & ~mTpAct;
    
    for i = find(q(:,1))' % go back to 1st increment
        cCkCrd{i} = cCkCrd{i}(mNoInc_grw(i,1):end,:);
        if mCkJun(i,2) ~= 0 % mod. intersection segment
            mCkJun(i,2) = mCkJun(i,2)-mNoInc_grw(i,1)+1;
        end
    end
    for i = find(q(:,2))' % go back to 1st increment
        cCkCrd{i} = cCkCrd{i}(1:end-mNoInc_grw(i,2)+1,:);
    end
    
    % restoring previous
    mTpRdi(q) = mTpRdi_ref(q);
    mTpAct(q) = true;
    mCkJun(q) = 0;
    
    % make unit growth
    mNoInc_grw(q) = 1;
    mNoInc_upd(q) = 1;
    
    % set new enrichment radii (not forgetting previous refinements)
    mTpRdi(mCk2Rf) = mTpRdi(mCk2Rf).*(0.5.^(nRefine_inc-mCkRfN(mCk2Rf)));
    
    % update n. times to refine
    mCkRfN(mCk2Rf) = nRefine_inc;
    
    % assert tips to be refined
    mCkRfn(mCk2Rf) = true;
    
    for i = find(mCk2Rf(:,1))'; d=cCkCrd{i}(1,:)-cCkCrd{i}(2,:);
        cCkCrd{i}(1,:) = cCkCrd{i}(2,:)+(f_inc*mTpRdi(i,1)/norm(d))*d;
    end
    for i = find(mCk2Rf(:,2))'; d=cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:);
        cCkCrd{i}(end,:) = cCkCrd{i}(end-1,:)+(f_inc*mTpRdi(i,2)/norm(d))*d;
    end
    
    % mNoInc_grw(mCk2Rf) = 1; % (by default)
    % mNoInc_upd(mCk2Rf) = 1; % (by default)
    
    % n.b.: "mNoInc_upd(mCk2Rf) = 1" is symbolic;
    % mNoInc_upd(mCk2Rf) = 0 after remeshing is done
    
end
end