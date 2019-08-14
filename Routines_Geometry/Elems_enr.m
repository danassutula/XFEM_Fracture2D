function [...
    vHvStd_elm,mHvStd_sgm,cHvBln_elm,cHvBln_sgm,cHvBln_rmp,cBrStd_eTp,...
    cBrStd_eSp,cBrStd_eFl,cBrBln_eSp,cBrBln_eFl,cBrBln_rSp,cBrBln_rFl...
    ] = Elems_enr(mNdCrd,mElNod,mCkCrd,vCkJun,r_tip,h0)

%--------------------------------------------------------------------------
% Initialise
%--------------------------------------------------------------------------

vHvStd_elm = [];
mHvStd_sgm = [];

cHvBln_elm = cell(1,2);
cHvBln_sgm = cell(1,2);
cHvBln_rmp = cell(1,2);

cBrStd_eSp = cell(1,2);
cBrBln_eSp = cell(1,2);
cBrBln_rSp = cell(1,2);

cBrStd_eFl = cell(1,2);
cBrBln_eFl = cell(1,2);
cBrBln_rFl = cell(1,2);

cBrStd_eTp = cell(1,2);

if isempty(mCkCrd)
   return 
end

tol = 1e-14;

if length(r_tip) == 1 
    r_tip(2) = r_tip;
end

%--------------------------------------------------------------------------
% Find cracked elements
%--------------------------------------------------------------------------

nCkCrd = size(mCkCrd,1);
nCkSgm = nCkCrd - 1;
nCkElm = nCkSgm*100; % estimate

vElCrk(nCkElm,1) = 0;
vCkSgm(nCkElm,1) = 0;

nCkElm = 0;

j_crk = [0;1];
r_tol = h0*tol;

for i_sgm = 1:nCkSgm
    
    j_crk = j_crk+1;
    x_sgm = mCkCrd(j_crk,:);
    
    % LS_tip: positive towards segment, negative outwards
    [LS_lin,LS_tp1,LS_tp2] = LevelSet2Line(mNdCrd,x_sgm);
    % LS_lin=LS_lin(:)';LS_tp1=LS_tp1(:)';LS_tp2=LS_tp2(:)';
    
    % elements whose nodes are on either side of the segment
    elmin = find(any(LS_lin(mElNod) < r_tol,2) & ... % below
                 any(LS_lin(mElNod) >-r_tol,2));     % above
    
    % elements who have at least one node in the interior of the segment
    elmin = elmin(any(ismember(mElNod(elmin,:),find(LS_tp1>-r_tol)),2) & ...
                  any(ismember(mElNod(elmin,:),find(LS_tp2>-r_tol)),2));
             
    if ~isempty(elmin) % (need this)
        
        if length(elmin) > 1
            p_end = any(LS_tp1(mElNod(elmin,:)) < r_tol,2) | ...
                    any(LS_tp2(mElNod(elmin,:)) < r_tol,2);
        else % wrong shape; need a row not a column vector
            p_end = any(LS_tp1(mElNod(elmin,:)) < r_tol) | ...
                    any(LS_tp2(mElNod(elmin,:)) < r_tol);
        end
        
        elmin_spl = elmin(~p_end); % definitely split
        elmin_end = elmin( p_end); % could be anything
        
        n = length(elmin_end);
        k = false(n,1);
        
        for i = 1:n
            k(i) = GeoXElem(x_sgm,mNdCrd(mElNod(elmin_end(i),:),:))>0;
        end
        
        % elements that are split or cut
        elmin = [elmin_spl;elmin_end(k)];
        
    end
    
    n = length(elmin);
    j = (1:n)+nCkElm;
    
    vElCrk(j) = elmin;
    vCkSgm(j) = i_sgm;
    
    nCkElm = nCkElm+n;
    
end

j = 1:nCkElm;
vElCrk = vElCrk(j); 
vCkSgm = vCkSgm(j);

[~,i_unq] = unique(vElCrk,'first'); i_unq = sort(i_unq);
i_rpt = true(length(vElCrk),1); i_rpt(i_unq) = false; 

% repeated (kink) elements (can have rpt. el)
vElRpt = vElCrk(i_rpt);
vSgRpt = vCkSgm(i_rpt);

% repeated (kink) elements (have no rpt. el)
[~,i_rpt] = unique(vElRpt,'last'); i_rpt = sort(i_rpt);
vElRpt = vElRpt(i_rpt); vSgRpt = vSgRpt(i_rpt);

% crack elements
vElCrk = vElCrk(i_unq); % unique crk. el.
mCkSgm = vCkSgm(i_unq,[1,1]); % init. crk. sgm. matching vElCrk

for i = 1:length(vElRpt) % modify mCkSgm(:,2)
    mCkSgm(vElRpt(i) == vElCrk,2) = vSgRpt(i);
end

%--------------------------------------------------------------------------
% Sort elements: std. and bln.
%--------------------------------------------------------------------------

if ~vCkJun(1) && ~vCkJun(2)

    %%% LEFT HAND SIDE OF CRACK
    
    vTpCrd = mCkCrd(1,:);
    jElTip = mCkSgm(:,1)==1;
    
    % find crack tip element
    for i = find(jElTip)'
        if ~GeoElm_xTip(vTpCrd,mNdCrd(mElNod(vElCrk(i),:),:))
            jElTip(i) = false;
        end
    end
    
    vElTip = vElCrk(jElTip);
    jElStd = ~jElTip; % !!! !!! !!!
    
    jElBln = find(jElStd); % init. Hvi. bln. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElTip,:));
    
    jBlRmp = sum(mBlRmp,2)>0; % id. of bln. el.
    jElBln = jElBln(jBlRmp); % bln. el. LHS
    jElStd(jElBln) = false; % upd std. el.
    
    cHvBln_elm{1} =  vElCrk(jElBln);
    cHvBln_sgm{1} =  mCkSgm(jElBln,:);
    cHvBln_rmp{1} = ~mBlRmp(jBlRmp,:);
    
    if r_tip(1)>0 && ~isempty(vElTip)
        
        [vElStd,vElBln,mBlRmp] = ...
            Elems_tip(mNdCrd,mElNod,vTpCrd,r_tip(1));
        
        if length(vElTip)==1
            jElTip = vElStd==vElTip; else
            jElTip = ismember(vElStd,vElTip);
        end
        
        cBrStd_eTp{1} = vElStd(jElTip);
        vElStd = vElStd(~jElTip);
        
        q = ismember(vElStd,vElCrk);
        cBrStd_eSp{1} = vElStd(q);
        
        q = ~q;
        cBrStd_eFl{1} = vElStd(q);
        
        q = ismember(vElBln,vElCrk);
        cBrBln_eSp{1} = vElBln(q);
        cBrBln_rSp{1} = mBlRmp(q,:);
        
        q = ~q;
        cBrBln_eFl{1} = vElBln(q);
        cBrBln_rFl{1} = mBlRmp(q,:);
        
    end
    
    %%% RIGHT HAND SIDE OF CRACK
    
    vTpCrd = mCkCrd(end,:);
    jElTip = mCkSgm(:,2)==size(mCkCrd,1)-1;
    
    % find crack tip element
    for i = find(jElTip)'
        if ~GeoElm_xTip(vTpCrd,mNdCrd(mElNod(vElCrk(i),:),:))
            jElTip(i) = false;
        end
    end
    
    vElTip = vElCrk(jElTip);
    jElStd = ~jElTip & jElStd; % !!! !!! !!!
    
    jElBln = find(jElStd); % init. Hvi. bln. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElTip,:));
    
    jBlRmp = sum(mBlRmp,2)>0; % id. of bln. el.
    jElBln = jElBln(jBlRmp); % bln. el. RHS
    jElStd(jElBln) = false; % upd std. el.
    
    cHvBln_elm{2} =  vElCrk(jElBln);
    cHvBln_sgm{2} =  mCkSgm(jElBln,:);
    cHvBln_rmp{2} = ~mBlRmp(jBlRmp,:); 
    
    if r_tip(2)>0 && ~isempty(vElTip)
        
        [vElStd,vElBln,mBlRmp] = ...
            Elems_tip(mNdCrd,mElNod,vTpCrd,r_tip(2));
        
        if length(vElTip)==1
            jElTip = vElStd==vElTip; else
            jElTip = ismember(vElStd,vElTip);
        end
        
        cBrStd_eTp{2} = vElStd(jElTip);
        vElStd = vElStd(~jElTip);
        
        q = ismember(vElStd,vElCrk);
        cBrStd_eSp{2} = vElStd(q);
        
        q = ~q;
        cBrStd_eFl{2} = vElStd(q);
        
        q = ismember(vElBln,vElCrk);
        cBrBln_eSp{2} = vElBln(q);
        cBrBln_rSp{2} = mBlRmp(q,:);
        
        q = ~q;
        cBrBln_eFl{2} = vElBln(q);
        cBrBln_rFl{2} = mBlRmp(q,:);
        
    end
    
elseif ~vCkJun(1) && vCkJun(2)  
    
    %%% LEFT HAND SIDE OF CRACK
    
    vTpCrd = mCkCrd(1,:);
    jElTip = mCkSgm(:,1)==1;
    
    % find crack tip element
    for i = find(jElTip)'
        if ~GeoElm_xTip(vTpCrd,mNdCrd(mElNod(vElCrk(i),:),:))
            jElTip(i) = false;
        end
    end
    
    vElTip = vElCrk(jElTip);
    jElStd = ~jElTip & mCkSgm(:,1)<vCkJun(2);
    
    jElBln = find(jElStd); % init. Hvi. bln. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElTip,:));
    
    jBlRmp = sum(mBlRmp,2)>0; % id. of bln. el.
    jElBln = jElBln(jBlRmp); % bln. el. LHS
    jElStd(jElBln) = false; % upd std. el.
    
    cHvBln_elm{1} =  vElCrk(jElBln);
    cHvBln_sgm{1} =  mCkSgm(jElBln,:);
    cHvBln_rmp{1} = ~mBlRmp(jBlRmp,:);
    
    if r_tip(1)>0 && ~isempty(vElTip)
        
        [vElStd,vElBln,mBlRmp] = ...
            Elems_tip(mNdCrd,mElNod,vTpCrd,r_tip(1));
        
        if length(vElTip)==1
            jElTip = vElStd==vElTip; else
            jElTip = ismember(vElStd,vElTip);
        end
        
        cBrStd_eTp{1} = vElStd(jElTip);
        vElStd = vElStd(~jElTip);
        
        q = ismember(vElStd,vElCrk);
        cBrStd_eSp{1} = vElStd(q);
        
        q = ~q;
        cBrStd_eFl{1} = vElStd(q);
        
        q = ismember(vElBln,vElCrk);
        cBrBln_eSp{1} = vElBln(q);
        cBrBln_rSp{1} = mBlRmp(q,:);
        
        q = ~q;
        cBrBln_eFl{1} = vElBln(q);
        cBrBln_rFl{1} = mBlRmp(q,:);
        
    end
    
    %%% RIGHT HAND SIDE OF CRACK
    
    jElBln = find(mCkSgm(:,1)>=vCkJun(2));
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),...
        mElNod(vElCrk(jElStd),:));
    
    jBlRmp = sum(mBlRmp,2)>0;
    jElBln = jElBln(jBlRmp); % bln. el. RHS
    % jElStd(jElBln) = false; % no need to upd.
    
    cHvBln_elm{2} = vElCrk(jElBln);
    cHvBln_sgm{2} = mCkSgm(jElBln,:);
    cHvBln_rmp{2} = mBlRmp(jBlRmp,:);

elseif vCkJun(1) && ~vCkJun(2)
    
    %%% RIGHT HAND SIDE OF CRACK
    
    vTpCrd = mCkCrd(end,:);
    jElTip = mCkSgm(:,2)==size(mCkCrd,1)-1;
    
    % find crack tip element
    for i = find(jElTip)'
        if ~GeoElm_xTip(vTpCrd,mNdCrd(mElNod(vElCrk(i),:),:))
            jElTip(i) = false;
        end
    end
    
    vElTip = vElCrk(jElTip);
    jElStd = ~jElTip & mCkSgm(:,2)>vCkJun(1);
    
    jElBln = find(jElStd); % init. Hvi. bln. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElTip,:));
    
    jBlRmp = sum(mBlRmp,2)>0; % id. of bln. el.
    jElBln = jElBln(jBlRmp); % bln. el. RHS
    jElStd(jElBln) = false; % upd std. el.
    
    cHvBln_elm{2} =  vElCrk(jElBln);
    cHvBln_sgm{2} =  mCkSgm(jElBln,:);
    cHvBln_rmp{2} = ~mBlRmp(jBlRmp,:); 
    
    if r_tip(2)>0 && ~isempty(vElTip)
        
        [vElStd,vElBln,mBlRmp] = ...
            Elems_tip(mNdCrd,mElNod,vTpCrd,r_tip(2));
        
        if length(vElTip)==1
            jElTip = vElStd==vElTip; else
            jElTip = ismember(vElStd,vElTip);
        end
        
        cBrStd_eTp{2} = vElStd(jElTip);
        vElStd = vElStd(~jElTip);
        
        q = ismember(vElStd,vElCrk);
        cBrStd_eSp{2} = vElStd(q);
        
        q = ~q;
        cBrStd_eFl{2} = vElStd(q);
        
        q = ismember(vElBln,vElCrk);
        cBrBln_eSp{2} = vElBln(q);
        cBrBln_rSp{2} = mBlRmp(q,:);
        
        q = ~q;
        cBrBln_eFl{2} = vElBln(q);
        cBrBln_rFl{2} = mBlRmp(q,:);
        
    end
    
    %%% LEFT HAND SIDE OF CRACK
    
    jElBln = find(mCkSgm(:,2)<=vCkJun(1));
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:), ...
        mElNod(vElCrk(jElStd),:));
    
    jBlRmp = sum(mBlRmp,2)>0;
    jElBln = jElBln(jBlRmp); % bln. el. LHS
    % jElStd(jElBln) = false; % no need to upd.
    
    cHvBln_elm{1} = vElCrk(jElBln);
    cHvBln_sgm{1} = mCkSgm(jElBln,:);
    cHvBln_rmp{1} = mBlRmp(jBlRmp,:);
    
else % vCkJun(1) && vCkJun(2) % two crack junctions (i.e. no crack tip)
    
    jElStd = mCkSgm(:,2)>vCkJun(1) ...
        & mCkSgm(:,1)<vCkJun(2); % with vtx. el.
    
    %%% LEFT HAND SIDE OF CRACK
    
    jElBln = find(mCkSgm(:,2)<=vCkJun(1) & ~jElStd); % bln. el. & no vtx. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElCrk(jElStd),:));
    
    jBlRmp = sum(mBlRmp,2)>0;
    jElBln = jElBln(jBlRmp);
    
    cHvBln_elm{1} = vElCrk(jElBln);
    cHvBln_sgm{1} = mCkSgm(jElBln,:);
    cHvBln_rmp{1} = mBlRmp(jBlRmp,:);
    
    %%% RIGHT HAND SIDE OF CRACK
    
    jElBln = find(mCkSgm(:,1)>=vCkJun(2) & ~jElStd); % bln. el. & no vtx. el.
    mBlRmp = ismember(mElNod(vElCrk(jElBln),:),mElNod(vElCrk(jElStd),:));
    
    jBlRmp = sum(mBlRmp,2)>0;
    jElBln = jElBln(jBlRmp);
    
    cHvBln_elm{2} = vElCrk(jElBln);
    cHvBln_sgm{2} = mCkSgm(jElBln,:);
    cHvBln_rmp{2} = mBlRmp(jBlRmp,:);
    
end
 
vHvStd_elm = vElCrk(jElStd);
mHvStd_sgm = mCkSgm(jElStd,:);

end