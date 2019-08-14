
%--------------------------------------------------------------------------
% Save current material state
%--------------------------------------------------------------------------

fprintf('\tbegin backup ...\n')

save_cBrStd_eTp = cBrStd_eTp;
save_cBrStd_eSp = cBrStd_eSp;
save_cBrStd_eFl = cBrStd_eFl;
save_cBrBln_eSp = cBrBln_eSp;
save_cBrBln_eFl = cBrBln_eFl;
save_cBrBln_rSp = cBrBln_rSp;
save_cBrBln_rFl = cBrBln_rFl;


save_cHvStd_elm = cHvStd_elm;
save_cHvStd_sgm = cHvStd_sgm;
save_cHvBln_elm = cHvBln_elm;
save_cHvBln_sgm = cHvBln_sgm;
save_cHvBln_rmp = cHvBln_rmp;


save_cNdEnr_brn = cNdEnr_brn;
save_cNdEnr_hvi = cNdEnr_hvi;
save_cNdBln_brn = cNdBln_brn;
save_cNdBln_hvi = cNdBln_hvi;


save_vElEnr = vElEnr;
save_cLNodE = cLNodE; 
save_cLEnDt = cLEnDt; 
save_cXsElm = cXsElm;


save_vElPhz = vElPhz;
save_mPrLod = mPrLod;


save_cGsEnr_omgShS = cGsEnr_omgShS;
save_cGsEnr_omgDvS = cGsEnr_omgDvS;
save_cGsEnr_omgShE = cGsEnr_omgShE;
save_cGsEnr_omgDvE = cGsEnr_omgDvE;
save_cGsEnr_omgWgt = cGsEnr_omgWgt;


save_nNdEnr = nNdEnr;
save_nGlDof = nGlDof;
save_nGDofE = nGDofE;
save_mNDofE = mNDofE;


if with_RfnXrs || with_RfnInc
   
    save_nNdStd = nNdStd;
    save_mNdCrd = mNdCrd;
    
    save_nElems = nElems;
    save_mLNodS = mLNodS;
    
    save_cElRfn = cElRfn;
    save_cElRf0 = cElRf0;
    
    save_cFxNod = cFxNod;
    save_cLdNod = cLdNod;
    save_cBnNod = cBnNod;
    
    % save_cBCNod = cBCNod;
    % save_vNdBnd = vNdBnd;
    
    save_vFxDof_all = vFxDof_all;
    save_vFxDof_nnz = vFxDof_nnz;
    
    % save all of Kg, Fg
    save_nGDofS = nGDofS;
    save_mGlStf = mGlStf;
    save_vGlFrc = vGlFrc;
    
else
    
    % only save Kse & Kee since symmetric
    save_mGlStf = mGlStf(:,nGDofS+1:end);
    save_vGlFrc = vGlFrc(nGDofS+1:end,1);
    
end

fprintf('\tbackup complete.\n\n')
