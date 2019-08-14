function [nnenr,elnod_enr] = Topo_omgEnr_add(nnenr,elnod_enr,elnew,elold,idcrk,idenr)

global mLNodS
global cLNodE
global cLEnDt

nnode_std = size(mLNodS,2);
nnode_enr = size(elnod_enr,2);
noenr = nnode_enr/nnode_std;

elnod_new = mLNodS(elnew,:);
elnod_old = mLNodS(elold,:);

% get member nodes (for glueing new el.)
[idmbr,idref] = ismember(elnod_new(:),elnod_old(:));

idmbr =  find(idmbr); % indices of member nodes (replace)
idref = idref(idmbr); % indices to member nodes (reuse)

q_rpt = ones(noenr,1); % (a trick for indexing)

idnod_new = ceil(idmbr/length(elnew));         % id of member nodes
idelm_new = idmbr-(idnod_new-1)*length(elnew); % id of member elements
idnod_new = idnod_new(:,q_rpt);                % (resize for mult. enr.)

[idref,~,q_m2r] = unique(idref); % q_m2r : i_mbr -> i_ref

idnod_old = ceil(idref/length(elold));         % id of member nodes
idelm_old = idref-(idnod_old-1)*length(elold); % id of member elements
idnod_old = idnod_old(:,q_rpt);                % (resize for mult. enr.)

for i = 2:noenr % adding periods if multiple enr.
    idnod_new(:,i) = idnod_new(:,i-1) + nnode_std;
    idnod_old(:,i) = idnod_old(:,i-1) + nnode_std;
end

n = numel(elnod_new);
q = idmbr(:,q_rpt);

% get enriched (blending) nodes to remove
for i = 2:noenr; q(:,i) = q(:,i-1)+n; end
ndbln = unique(elnod_enr(q(:))); % sorted !!!
nndel = length(ndbln); % n. nd. to delete

% shift enr. nd.
for i = nndel:-1:1
    q = elnod_enr(:) > ndbln(i);
    elnod_enr(q) = elnod_enr(q)-1;
end

% enr. (bln.) nd. to reuse
ndbln = idnod_old; % init.

% short-list enr. nd. to reuse in new enr. el. topo.
for i = 1:length(idref); k = elold(idelm_old(i));
    
    % determine node position within enr. data
    q = cLEnDt{k}(1,:) == idcrk & cLEnDt{k}(2,:) == idenr;
    if isempty(q); error('could not find enrichment.'); end
    n = sum(cLEnDt{k}(3,1:find(q,1,'first')),2);
    q = n+1-nnode_enr:n;
    
    % store enr. nodes to be reused
    ndbln(i,:) = cLNodE{k}(q(idnod_old(i,:)));
    
end

for i = 1:length(idmbr)
    elnod_enr(idelm_new(i),idnod_new(i,:)) = ndbln(q_m2r(i),:);
end

% update n. enr. nodes
nnenr = nnenr - nndel;

% if nargout == 1
%     for j = 1:length(elnew) % store enr. nodes (push back)
%         cLNodE{elnew(j)}(end+1:end+nnode_enr) = elnod_enr(j,:);
%     end
% end

end