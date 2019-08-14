function [ndenr,enunq] = Topo_omgEnr_rmv(elenr)

% ndenr - enriched nodes that have been removed, size = [N,1], type = logical
% enunq - unique enrichments that existed in elements 'elenr', size = [N,2] 

if isempty(elenr)
    ndenr = [];
    enunq = []; % (may need empty but with 2 col.)
    return
end

global mLNodS mNDofE cLNodE cLEnDt cXsElm       ...
    nEnFun_brn cBrStd_eTp cBrStd_eSp cBrStd_eFl ...
    cBrBln_eSp cBrBln_eFl cBrBln_rSp cBrBln_rFl ...
    cHvStd_elm cHvBln_elm cHvStd_sgm cHvBln_sgm ...
    cHvBln_rmp cNdEnr_brn cNdEnr_hvi cNdBln_hvi

nNdEn8 = size(mNDofE,2);
nLNodS = size(mLNodS,2);
nLNodE_brn = nEnFun_brn*nLNodS;

enunq = cat(2,cLEnDt{elenr});
enunq = enunq(1:2,:)'; % get types of enr.
enunq = unique(enunq,'rows'); % unq. enr.

% enr. nd. to remove
ndenr = false(nNdEn8,1);
nlenr = length(elenr);

for i = 1:size(enunq,1) % process all enrichments
    
    ndths = false(nNdEn8,1);
    elths = false(nlenr,1);
    
    for k = 1:nlenr % get el. in connection to this (i) enrichment
        elths(k) = any(cLEnDt{elenr(k)}(1,:) == enunq(i,1) & ...
            cLEnDt{elenr(k)}(2,:) == enunq(i,2));
    end
    
    % el. with this enr.
    elths = elenr(elths);
    
    % this crack
    k = enunq(i,1);
    switch enunq(i,2)
        case 1 % (brn. tip #1)
            
            % all elements that are affected by this enrichment
            elth0 = [cBrStd_eTp{k,1};cBrStd_eSp{k,1};...
                cBrStd_eFl{k,1};cBrBln_eSp{k,1};cBrBln_eFl{k,1}];
            
            % elements outside zone
            elth0 = setdiff(elth0,elths);
            
            % bool to enr. nd. to del. (n.b.: mult. enr. fun.)
            nombr = ~ismember(mLNodS(elths,:),mLNodS(elth0,:));
            q=(1:nLNodS)'; nombr=nombr(:,q(:,ones(1,nEnFun_brn)));
            
            for jj = 1:length(elths); j = elths(jj);
                
                % determine position within enr. data of this enr.
                q = cLEnDt{j}(1,:) == k & cLEnDt{j}(2,:) == 1;
                n = sum(cLEnDt{j}(3,1:find(q,1,'first')),2);
                
                q = n+1-nLNodE_brn:n; % all el. enr. nodes of this enr.
                ndths(cLNodE{j}(q(nombr(jj,:)))) = true;
                
            end
            
            p = find(ndths);
            
            % update array of enr. nodes for this (Brn.) enr.
            cNdEnr_brn{k,1} = setdiff(cNdEnr_brn{k,1},p);
            % cNdBln_brn{k,1} = setdiff(cNdBln_brn{k,1},p);
            
            cBrStd_eTp{k,1} = setdiff(cBrStd_eTp{k,1},elenr);
            cBrStd_eSp{k,1} = setdiff(cBrStd_eSp{k,1},elenr);
            cBrStd_eFl{k,1} = setdiff(cBrStd_eFl{k,1},elenr);
            
            q = ~ismember(cBrBln_eSp{k,1},elenr);
            cBrBln_eSp{k,1} = cBrBln_eSp{k,1}(q);
            cBrBln_rSp{k,1} = cBrBln_rSp{k,1}(q,:);
            
            q = ~ismember(cBrBln_eFl{k,1},elenr);
            cBrBln_eFl{k,1} = cBrBln_eFl{k,1}(q);
            cBrBln_rFl{k,1} = cBrBln_rFl{k,1}(q,:);
            
        case 2 % (brn. tip #2)
            
            % all elements that are affected by this enrichment
            elth0 = [cBrStd_eTp{k,2};cBrStd_eSp{k,2};...
                cBrStd_eFl{k,2};cBrBln_eSp{k,2};cBrBln_eFl{k,2}];
            
            % elements outside zone
            elth0 = setdiff(elth0,elths);
            
            % bool to enr. nd. to del. (n.b.: mult. enr. fun.)
            nombr = ~ismember(mLNodS(elths,:),mLNodS(elth0,:));
            q=(1:nLNodS)'; nombr=nombr(:,q(:,ones(1,nEnFun_brn)));
            
            for jj = 1:length(elths); j = elths(jj);
                
                % determine position within enr. data of this enr.
                q = cLEnDt{j}(1,:) == k & cLEnDt{j}(2,:) == 2;
                n = sum(cLEnDt{j}(3,1:find(q,1,'first')),2);
                
                q = n+1-nLNodE_brn:n; % all el. enr. nodes of this enr.
                ndths(cLNodE{j}(q(nombr(jj,:)))) = true;
                
            end
            
            p = find(ndths);
            
            % update array of enr. nodes for this (Brn.) enr.
            cNdEnr_brn{k,2} = setdiff(cNdEnr_brn{k,2},p);
            % cNdBln_brn{k,2} = setdiff(cNdBln_brn{k,2},p));
            
            cBrStd_eTp{k,2} = setdiff(cBrStd_eTp{k,2},elenr);
            cBrStd_eSp{k,2} = setdiff(cBrStd_eSp{k,2},elenr);
            cBrStd_eFl{k,2} = setdiff(cBrStd_eFl{k,2},elenr);
            
            q = ~ismember(cBrBln_eSp{k,2},elenr);
            cBrBln_eSp{k,2} = cBrBln_eSp{k,2}(q);
            cBrBln_rSp{k,2} = cBrBln_rSp{k,2}(q,:);
            
            q = ~ismember(cBrBln_eFl{k,2},elenr);
            cBrBln_eFl{k,2} = cBrBln_eFl{k,2}(q);
            cBrBln_rFl{k,2} = cBrBln_rFl{k,2}(q,:);
            
        otherwise % == 0 (Hvi.)
            
            % all elements that are affected by this enrichment
            elth0 = [cHvBln_elm{k,1};cHvStd_elm{k};cHvBln_elm{k,2}];
            elth0 = setdiff(elth0,elths); % elements outside zone
            
            % bool to enr. nd. to del. (n.b.: single enr. fun.)
            nombr = ~ismember(mLNodS(elths,:),mLNodS(elth0,:));
            % q=1:nLNodS; nombr=nombr(:,q(:,ones(1,nEnFun_hvi)));
            
            for jj = 1:length(elths); j = elths(jj);
                
                % determine position within enr. data of this enr.
                q = cLEnDt{j}(1,:) == k & cLEnDt{j}(2,:) == 0;
                n = sum(cLEnDt{j}(3,1:find(q,1,'first')),2);
                
                q = n+1-nLNodS:n; % all el. enr. nodes of this enr.
                ndths(cLNodE{j}(q(nombr(jj,:)))) = true;
                
            end
            
            p = find(ndths);
            
            % update array of enr. nodes for this (Hvi.) enr.
            cNdEnr_hvi{k} = setdiff(cNdEnr_hvi{k},p);
            
            % update array of bln. nodes (do both ends to be safe)
            cNdBln_hvi{k,1} = setdiff(cNdBln_hvi{k,1},p);
            cNdBln_hvi{k,2} = setdiff(cNdBln_hvi{k,2},p);
            
            q = ~ismember(cHvStd_elm{k},elenr);
            cHvStd_elm{k} = cHvStd_elm{k}(q);
            cHvStd_sgm{k} = cHvStd_sgm{k}(q,:);
            
            q = ~ismember(cHvBln_elm{k,1},elenr);
            cHvBln_elm{k,1} = cHvBln_elm{k,1}(q);
            cHvBln_sgm{k,1} = cHvBln_sgm{k,1}(q,:);
            cHvBln_rmp{k,1} = cHvBln_rmp{k,1}(q,:);
            
            q = ~ismember(cHvBln_elm{k,2},elenr);
            cHvBln_elm{k,2} = cHvBln_elm{k,2}(q);
            cHvBln_sgm{k,2} = cHvBln_sgm{k,2}(q,:);
            cHvBln_rmp{k,2} = cHvBln_rmp{k,2}(q,:);
            
    end
    
    % upd. enr. nd. to rmv.
    ndenr(ndths) = true;
    
end

for i = elenr(:)'
    cLNodE{i} = [];
    cLEnDt{i} = [];
    cXsElm{i} = [];
end
