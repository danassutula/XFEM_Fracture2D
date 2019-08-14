function [ndEnr,lNods_enr] = Topo_omgEnr(ndEnr,noEnr,lNods_ref)
%==========================================================================
%{
    OUTPUT:
    ndEnr(1)       = updated ref. enr. node
    lNods_enr(:,:) = new element topology (elements nodes are in order)

    INPUT:
    ndEnr(1)       = ref. enr. node
    noEnr(1)       = no. enr functions
    lNods_ref(:,:) = element topology (sub)
%}
%==========================================================================
%                               BULLETPROOF
%==========================================================================

[n_elm,n_nod] = size(lNods_ref);
lNods_enr = zeros(n_elm,n_nod*noEnr);
[ndRef,lNods_ref] = Topo_omgRef(0,lNods_ref);

j_nod = 1:n_nod;
for i = 1:noEnr
    
    lNods_enr(:,j_nod) = lNods_ref + ndEnr;
    
    ndEnr = ndEnr + ndRef;
    j_nod = j_nod + n_nod;
    
end

%==========================================================================

end