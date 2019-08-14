
%--------------------------------------------------------------------------
% Restore previous material state
%--------------------------------------------------------------------------

fprintf('\tbegin restore ...\n')

cBrStd_eTp = save_cBrStd_eTp;
cBrStd_eSp = save_cBrStd_eSp;
cBrStd_eFl = save_cBrStd_eFl;
cBrBln_eSp = save_cBrBln_eSp;
cBrBln_eFl = save_cBrBln_eFl;
cBrBln_rSp = save_cBrBln_rSp;
cBrBln_rFl = save_cBrBln_rFl;


cHvStd_elm = save_cHvStd_elm;
cHvStd_sgm = save_cHvStd_sgm;
cHvBln_elm = save_cHvBln_elm;
cHvBln_sgm = save_cHvBln_sgm;
cHvBln_rmp = save_cHvBln_rmp;


cNdEnr_brn = save_cNdEnr_brn;
cNdEnr_hvi = save_cNdEnr_hvi;
cNdBln_brn = save_cNdBln_brn;
cNdBln_hvi = save_cNdBln_hvi;


vElEnr = save_vElEnr;
cLNodE = save_cLNodE; 
cLEnDt = save_cLEnDt; 
cXsElm = save_cXsElm;


vElPhz = save_vElPhz;
mPrLod = save_mPrLod;


cGsEnr_omgShS = save_cGsEnr_omgShS;
cGsEnr_omgDvS = save_cGsEnr_omgDvS;
cGsEnr_omgShE = save_cGsEnr_omgShE;
cGsEnr_omgDvE = save_cGsEnr_omgDvE;
cGsEnr_omgWgt = save_cGsEnr_omgWgt;


nNdEnr = save_nNdEnr;
nGlDof = save_nGlDof;
nGDofE = save_nGDofE;
mNDofE = save_mNDofE;


if with_RfnXrs || with_RfnInc
    
    nNdStd = save_nNdStd;
    mNdCrd = save_mNdCrd;
    
    nElems = save_nElems;
    mLNodS = save_mLNodS;
    
    cElRfn = save_cElRfn;
    cElRf0 = save_cElRf0;
    
    cFxNod = save_cFxNod;
    cLdNod = save_cLdNod;
    cBnNod = save_cBnNod;
    
    % cBCNod = save_cBCNod;
    % vNdBnd = save_vNdBnd;
    
    vNdBnd = cat(1,cBnNod{:});
    vNdBnd = vNdBnd(:,1);
    
    vFxDof_all = save_vFxDof_all;
    vFxDof_nnz = save_vFxDof_nnz;
    
    nGDofS = save_nGDofS;
    mGlStf = save_mGlStf;
    vGlFrc = save_vGlFrc;
    
else
    
    if nGlDof < length(vGlFrc)
        mGlStf(:,nGlDof+1:end) = [];
        mGlStf(nGlDof+1:end,:) = [];
        vGlFrc(nGlDof+1:end,:) = [];
    else
        mGlStf(nGlDof,nGlDof) = 0;
        vGlFrc(nGlDof,     1) = 0;
    end
    
    mGlStf(:,nGDofS+1:end) = save_mGlStf;
    mGlStf(nGDofS+1:end,:) = save_mGlStf'; % (overlapping, but faster!)
    vGlFrc(nGDofS+1:end,1) = save_vGlFrc;

end
    
fprintf('\trestore complete.\n\n')
