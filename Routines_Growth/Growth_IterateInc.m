
%==========================================================================
% Iterations: global energy minimization
%==========================================================================


%--------------------------------------------------------------------------
% Assess the energy release rate at each (incremented) crack tip
%--------------------------------------------------------------------------

R_sif = mTpRdi*f_sif; % 'mTpRdi' is post-intersection/refinement
R_sif(~mCk2Up_itr)=0; % 'mCk2Up_itr' is pre-intersection

Ji = IntJ(cCkCrd,mNDspS,mNDspE, ...
    mPrLod,vElPhz,cDMatx,R_sif);

if wCkLod
    Ji = Ji + IntJ_CrkLodHvi(cCkCrd,mNDspE, ...
        nGsHvi_gam,mGsHvi_gamShp,vGsHvi_gamWgt,mCkLod,R_sif);
    Ji = Ji + IntJ_CrkLodBrn(cCkCrd,mNDspE, ...
        nGsBrn_gam,mGsBrn_gamShp,vGsBrn_gamWgt,mCkLod,R_sif);
end

if wBdLod
    warning('body force ignored in J-integral.')
end

%--------------------------------------------------------------------------
% Identify crack tips to discard
%--------------------------------------------------------------------------

if exist('RmvCritr_iterInc','var')
    q = mTpAct & mCk2Up_itr & RmvCritr_iterInc(Ji);
else
    q = mTpAct & mCk2Up_itr & Ji<sum(Ji(:))/nnz(Ji) * 0.98;
end

%--------------------------------------------------------------------------
% Discard crack tips below threshold
%--------------------------------------------------------------------------

if any(q(:)) % reestablish previous crack state
    
    for i = find(q(:,1))'; cCkCrd{i}=cCkCrd{i}(2:end,:);
        if mCkJun(i,2); mCkJun(i,2) = mCkJun(i,2)-1; end
    end
    
    for i = find(q(:,2))'
        cCkCrd{i} = cCkCrd{i}(1:end-1,:);
    end
    
    if with_RfnXrs || with_AdpEnr
        mTpRdi(q) = mTpRdi_ref(q);
    end
    
    if with_RfnXrs % || with_RfnInc
        mCkRfn(q) = mCkRfn_ref(q);
        mCkRfN(q) = mCkRfN_ref(q);
        mIsInc(q) = mIsInc_ref(q);
    end
    
    mCk2Up_itr(q) = 0;
    mCk2Up = mCk2Up_itr;
    
    vTp2Up = find(mCk2Up(:));
    nTp2Up = length(vTp2Up);
    
    mNoInc_grw = real(mCk2Up);
    mNoInc_upd = mNoInc_grw;
    
    if with_GLwDir
        fGlLaw_dir = 1; % (unconditional)
    end
    
else
    
    fIter_inc = 0;
    mCk2Up(:) = 0;
    
    vTp2Up = [];
    nTp2Up = 0;
    
    mNoInc_grw(:) = 0;
    mNoInc_upd(:) = 0;
    
end

%--------------------------------------------------------------------------
