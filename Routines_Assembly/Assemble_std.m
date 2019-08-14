
%==========================================================================
% Assemble global (std.)
%==========================================================================


%--------------------------------------------------------------------------
% Assemble global stiffness matrix
%--------------------------------------------------------------------------

nGDofS = nNdStd*nDimes;

if with_RdoStd || ~exist('mGStfS','var')
    mGStfS = MtxStiffStd(nGDofS,mNdCrd,mLNodS,...
        vElPhz,cDMatx,mGsStd_omgDrv,vGsStd_omgWgt);
else
    if size(mGStfS,1) ~= nGDofS
        warning('wrong size of stiffness matrix; re-doing')
        mGStfS = MtxStiffStd(nGDofS,mNdCrd,mLNodS,...
            vElPhz,cDMatx,mGsStd_omgDrv,vGsStd_omgWgt);
    end
end

%--------------------------------------------------------------------------
% Singular nodes
%--------------------------------------------------------------------------

vSgDof_std = (1:nGDofS)'; % initialise
vSgDof_std = vSgDof_std(diag(mGStfS(vSgDof_std,vSgDof_std)) == 0);
nSgDof_std = length(vSgDof_std);

if ~isempty(vSgDof_std)
   sprintf('Std. singular DOFs detected: n = %d\n',nSgDof_std)
end

%--------------------------------------------------------------------------
% Initialise global force vector
%--------------------------------------------------------------------------

vGFrcS = sparse([],[],[],nGDofS,1);

%--------------------------------------------------------------------------
% Boundary tractions
%--------------------------------------------------------------------------

for i = 1:nBnLod; bnlod = mBnLod(i,:);
    if any(bnlod)
        if size(cLdNod{i},2) > 1
            
            % 'line' load
            vGFrcS = vGFrcS + ForceStd_Bound(nGDofS,mNdCrd,cLdNod{i},...
                mGsStd_gamShp,mGsStd_gamDrv,vGsStd_gamWgt,bnlod);
            
        else
            
            % 'point' load
            k = cLdNod{i}*nDimes; j = k-nDimes+1:k; 
            vGFrcS(j) = vGFrcS(j) + bnlod(:);
            
        end
    end
end

%--------------------------------------------------------------------------
% Body tractions
%--------------------------------------------------------------------------

if wBdLod
    
    vGFrcS = vGFrcS + ForceStd_Body(nGDofS,mNdCrd,mLNodS,...
        mGsStd_omgShp,mGsStd_omgDrv,vGsStd_omgWgt,vBdLod);
    
end

%--------------------------------------------------------------------------
% Pre-stress/pre-strain/thermal
%--------------------------------------------------------------------------

if wPrLod
    
    vGFrcS = vGFrcS + ForceStd_Resid(nGDofS,mNdCrd,mLNodS,...
        mGsStd_omgDrv,vGsStd_omgWgt,mPrLod);
    
end

%--------------------------------------------------------------------------
% Clear memory
%--------------------------------------------------------------------------

clear bnlod

%--------------------------------------------------------------------------
