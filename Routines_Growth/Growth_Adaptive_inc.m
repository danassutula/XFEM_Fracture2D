
%--------------------------------------------------------------------------
% Adaptive fracture advance/remeshing based on increment kinking
%--------------------------------------------------------------------------

if with_GLwDir % (or if with_DirAvg == 1)
if dBetaMax_mshCors ~= 0 % (no coarsening allowed)
    
    % ADAPTIVE COARSENING
    
    % n.b. useless to test against abs(mTpBta)<dBetaMax_mshCors if
    % running pure energy minimization algo. (i.e. with_DirAvg == 0);
    % this is because the previously converged increment direction is
    % predominantly mode-I such that for the current tips K_II/K_I->0
    
    % initial guess of which crack tips to coarsen (mCk2Rf is local)
    mCk2Rf = mNoInc_grw>0 & mTpAct & mCkRfN>0; % & abs(mTpBta)<dBetaMax_mshCors;
    
    % to coarsen if kink of prev. inc. was small enough
    for i = find(mCk2Rf(:,1))'; x = cCkCrd{i}(2:4,:);
        mCk2Rf(i,1) = ~(abs(atan2(x(1,2)-x(2,2),x(1,1)-x(2,1)) - ...
            atan2(x(2,2)-x(3,2),x(2,1)-x(3,1)))>dBetaMax_mshCors);
    end
    
    for i = find(mCk2Rf(:,2))'; x = cCkCrd{i}(end-3:end-1,:);
        mCk2Rf(i,2) = ~(abs(atan2(x(3,2)-x(2,2),x(3,1)-x(2,1)) - ...
            atan2(x(2,2)-x(1,2),x(2,1)-x(1,1)))>dBetaMax_mshCors);
    end
    
    if any(mCk2Rf(:))
        
        % update tip positions: doubling the increment
        fprintf('\n\tadaptive remeshing: coarsening\n')
        
        for i = find(mCk2Rf(:,1))'; pass=1;
            for j = setdiff(1:nCrack,i)
                % if min. dist. is violated by doubling the inc. length
                if LevelSet2Poly_abs(cCkCrd{i}(1,:),cCkCrd{j}) < 3*mTpRdi(i,1)
                    pass = 0; mCk2Rf(i,1) = false; break; % flag not to coarsen
                end
            end
            if pass == 1; cCkCrd{i}(1,:) = ... % double increment length
                cCkCrd{i}(2,:) + 2*diff(cCkCrd{i}([2,1],:),1);
            end
        end
        
        for i = find(mCk2Rf(:,2))'; pass=1;
            for j = setdiff(1:nCrack,i)
                % if min. dist. is violated by doubling the inc. length
                if LevelSet2Poly_abs(cCkCrd{i}(end,:),cCkCrd{j}) < 3*mTpRdi(i,2)
                    pass = 0; mCk2Rf(i,2) = false; break; % flag not to coarsen
                end
            end
            if pass == 1; cCkCrd{i}(end,:) =  ... % double increment length
                cCkCrd{i}(end-1,:) + 2*diff(cCkCrd{i}([end-1,end],:),1);
            end
        end
        
        % update n. times to refine
        mCkRfN = mCkRfN-real(mCk2Rf);
        
        % set new enrichment radii: doubling
        mTpRdi = mTpRdi.*(2.^real(mCk2Rf));
        
        % mNoInc_grw(mCk2Rf) = 1; % (by default)
        % mNoInc_upd(mCk2Rf) = 1; % (by default)
        
        % n.b.: "mNoInc_upd(mCk2Rf) = 1" is symbolic;
        % mNoInc_upd(mCk2Rf) = 0 after remeshing is done
        
    end
end
else % with_GLwDir == 0
    
    % ADAPTIVE REFINEMENT
    
    % n.b. mIsInc ensures that a specific crack has already propagated;
    % only then refinement and crack increment halfing can be carried out;
    % otherwise, halfing of the increment would undermine the initial crack
    
    mCk2Rf = mNoInc_grw>0 & mCkRfN<nRefine_inc & ...
        abs(mTpBta)>dAngleMin_msh & mIsInc; % (mCk2Rf is local)
    
    if any(mCk2Rf(:))
        
        % update tip positions: halfing the increment
        fprintf('\n\tadaptive remeshing: refining\n')
        
        % pull back intersections
        q = mCk2Rf & ~mTpAct;
        
        for i = find(q(:,1))' % go back to 1st increment
            cCkCrd{i} = cCkCrd{i}(mNoInc_grw(i,1):end,:);
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
        
        % assert tips to be refined
        mCkRfn(mCk2Rf) = true;
        
        % update n. times to refine
        mCkRfN = mCkRfN+real(mCk2Rf);
        
        % set new enrichment radii: halfing
        mTpRdi = mTpRdi.*(0.5.^real(mCk2Rf));
        
        for i = find(mCk2Rf(:,1))'
            cCkCrd{i} = cCkCrd{i}(2:end,:); % go back 1 inc.
            cCkCrd{i}(1,:) = 0.5*sum(cCkCrd{i}([1,2],:),1);
        end
        for i = find(mCk2Rf(:,2))'
            cCkCrd{i} = cCkCrd{i}(1:end-1,:); % go back 1 inc.
            cCkCrd{i}(end,:) = 0.5*sum(cCkCrd{i}([end-1,end],:),1);
        end
        
        mNoInc_grw(mCk2Rf) = 0;
        % mNoInc_upd(mCk2Rf) = 1; % (by default)
        
        % n.b.: "mNoInc_upd(mCk2Rf) = 1" is symbolic;
        % mNoInc_upd(mCk2Rf) = 0 after remeshing is done
        
    end
    
    % ADAPTIVE COARSENING
    
    mCk2Rf = mNoInc_grw>0 & mCkRfN>0 & ...
        abs(mTpBta)<dBetaMax_mshCors & mTpAct; % (mCk2Rf is local)
    
    if any(mCk2Rf(:))
        
        % update tip positions: doubling the increment
        fprintf('\n\tadaptive remeshing: coarsening\n')
        
        for i = find(mCk2Rf(:,1))'; pass=1;
            for j = setdiff(1:nCrack,i)
                % check if min. dist. is violated by doubling the inc. length
                if LevelSet2Poly_abs(cCkCrd{i}(1,:),cCkCrd{j}) < 3*mTpRdi(i,1)
                    pass = 0; mCk2Rf(i,1) = false; break; % flag not to coarsen
                end
            end
            if pass == 1; cCkCrd{i}(1,:) = ... % double increment length
                cCkCrd{i}(2,:) + 2*diff(cCkCrd{i}([2,1],:),1);
            end
        end
        
        for i = find(mCk2Rf(:,2))'; pass=1;
            for j = setdiff(1:nCrack,i)
                % check if min. dist. is violated by doubling the inc. length
                if LevelSet2Poly_abs(cCkCrd{i}(end,:),cCkCrd{j}) < 3*mTpRdi(i,2)
                    pass = 0; mCk2Rf(i,2) = false; break; % flag not to coarsen
                end
            end
            if pass == 1; cCkCrd{i}(end,:) =  ... % double increment length
                cCkCrd{i}(end-1,:) + 2*diff(cCkCrd{i}([end-1,end],:),1);
            end
        end
        
        % update n. times to refine
        mCkRfN = mCkRfN-real(mCk2Rf);
        
        % set new enrichment radii: doubling
        mTpRdi = mTpRdi.*(2.^real(mCk2Rf));
        
        % mNoInc_grw(mCk2Rf) = 1; % (by default)
        % mNoInc_upd(mCk2Rf) = 1; % (by default)
        
        % n.b.: "mNoInc_upd(mCk2Rf) = 1" is symbolic;
        % mNoInc_upd(mCk2Rf) = 0 after remeshing is done
        
    end
end