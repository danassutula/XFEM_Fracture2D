function [cCkCrd,mTpAct,mTpRdi,mCkJun] = ...
    GeoCrk_Intersect0_mindist_all(cCkCrd,mTpAct,mTpRdi,mCkJun,f_tip,f_xrs)

% checks a particular crack (iCrack) against all cracks and intersects 
% the one (jCrack) that is closest and satisfies the min. dist. criterion

% min dist. to ck. vtx. as frac. of 'he'
f_vtx = 0.1;

% max no. bln. sgm. to add
nSgAdd = 2;

msg1 = 'complex crack topology; removing tip enrichment; i_crk = %d, i_tip = 1';
msg2 = 'complex crack topology; removing tip enrichment; i_crk = %d, i_tip = 2';

% rel. to 'f_tip'
f_vtx = f_vtx/f_tip;

nCrack = length(cCkCrd);
vCkAll = nCrack:-1:1;

for iCrack = vCkAll % intersecting crack
    
    if mTpAct(iCrack,1)
        
        %------------------------------------------------------------------
        % TIP #1
        %------------------------------------------------------------------
        
        jCkSgm = mCkJun(iCrack,:);
        nCkCrd = size(cCkCrd{iCrack},1);
        
        if ~jCkSgm(2)
            jCkSgm(2) = nCkCrd;
        end
        
        mCkCrd_ref = cCkCrd{iCrack}(1:jCkSgm(2),:);
        vTpCrd_ref = mCkCrd_ref(1,:);
        
        % tolerance for min. dist.
        tol_xrs = mTpRdi(iCrack,1)*f_xrs;
        tol_vtx = mTpRdi(iCrack,1)*f_vtx;
        
        % keep track of min-dist. intersections for iCrack
        mindist_r = zeros(nCrack,1); % radius to jCrack
        mindist_d = zeros(nCrack,2); % vector to jCrack
        mindist_i = zeros(nCrack,1); % segment of jCrack
        
        mindist_r(iCrack) = Inf;
        
        for jCrack = vCkAll(vCkAll ~= iCrack) % intersected crack (maybe)
            
            jCkSgm = mCkJun(jCrack,:);
            nCkCrd = size(cCkCrd{jCrack},1);
            
            if jCkSgm(1)
                jCkSgm(1) = jCkSgm(1)+1;
            else
                jCkSgm(1) = 1;
            end
            
            if ~jCkSgm(2)
                jCkSgm(2) = nCkCrd;
            end
            
            mCkCrd = cCkCrd{jCrack}(jCkSgm(1):jCkSgm(2),:);
            nCkCrd = jCkSgm(2)-jCkSgm(1)+1;
            
            % get distance to crack
            [r,d_inc,i_sgm] = LevelSet2Poly_abs(vTpCrd_ref,mCkCrd);
            
            mindist_r(jCrack,1) = r;
            mindist_d(jCrack,:) = d_inc;
            mindist_i(jCrack,1) = i_sgm;
            
        end
        
        [r,i] = min(mindist_r);
        d_inc = mindist_d(i,:);
        i_sgm = mindist_i(i,:);
        
        jCrack = i; % min-dist. crack
        jCkSgm = mCkJun(jCrack,:);
        nCkCrd = size(cCkCrd{jCrack},1);
        
        if jCkSgm(1)
            jCkSgm(1) = jCkSgm(1)+1;
        else
            jCkSgm(1) = 1;
        end
        
        if ~jCkSgm(2)
            jCkSgm(2) = nCkCrd;
        end
        
        mCkCrd = cCkCrd{jCrack}(jCkSgm(1):jCkSgm(2),:);
        nCkCrd = jCkSgm(2)-jCkSgm(1)+1;
        
        if r < tol_xrs
            
            % tip segment direction
            d_tip = vTpCrd_ref - mCkCrd_ref(2,:);
            
            if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1)) ...
                    && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0 % acute
                
                % A ----x--> B
                
                x_crs = vTpCrd_ref + d_inc;
                
                x_ptA = mCkCrd(i_sgm,:);
                x_ptB = mCkCrd(i_sgm+1,:);
                
                d_sgm = x_ptB - x_ptA;
                d_xsg = x_crs - x_ptA;
                
                s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                
                % i_sgm is index into mCkCrd, then
                % j_sgm will be index into cCkCrd{jCrack}
                j_sgm = jCkSgm(1)+i_sgm-1;
                
                if s_xsg < tol_vtx % xs is on A
                    
                    if i_sgm ~= 1 || ~mCkJun(jCrack,1)
                        if i_sgm == 1 % xs is on tip
                            
                            mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                            mTpAct(jCrack,1) = false;
                            mTpRdi(jCrack,1) = 0;
                            
                        else
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                % mCkAdd = mCkCrd(end:-1:i_sgm,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm,:);
                            else
                                % mCkAdd = mCkCrd(1:i_sgm,:);
                                mCkAdd = cCkCrd{jCrack}(1:j_sgm,:);
                            end
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(end-nSgAdd:end,:);
                        end
                        
                        mCkJun(iCrack,1) = size(mCkAdd,1)-1;
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                        
                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                                mCkJun(iCrack,2) + size(mCkAdd,1);
                        end
                        
                        mTpAct(iCrack,1) = false;
                        mTpRdi(iCrack,1) = 0;
                        
%                         break
                        
                    end
                    
                elseif s_sgm-s_xsg < tol_vtx % xs on B
                    
                    if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                        if i_sgm+1 == nCkCrd % xs is on tip
                            
                            mCkAdd = cCkCrd{jCrack};
                            mTpAct(jCrack,2) = false;
                            mTpRdi(jCrack,2) = 0;
                            
                        else
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                % mCkAdd = mCkCrd(end:-1:i_sgm+1,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm+1,:);
                            else
                                % mCkAdd = mCkCrd(1:i_sgm+1,:);
                                mCkAdd = cCkCrd{jCrack}(1:j_sgm+1,:);
                            end
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(end-nSgAdd:end,:);
                        end
                        
                        mCkJun(iCrack,1) = size(mCkAdd,1)-1;
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                        
                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                                mCkJun(iCrack,2) + size(mCkAdd,1);
                        end
                        
                        mTpAct(iCrack,1) = false;
                        mTpRdi(iCrack,1) = 0;
                        
%                         break
                        
                    end
                    
                else % xs is between A and B
                    
                    if nCkCrd > 2
                        if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                            % mCkAdd = [mCkCrd(end:-1:i_sgm+1,:);x_crs];
                            mCkAdd = [cCkCrd{jCrack}(end:-1:j_sgm+1,:);x_crs];
                        else
                            % mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                            mCkAdd = [cCkCrd{jCrack}(1:j_sgm,:);x_crs];
                        end
                    else
                        if s_sgm > 2*s_xsg
                            mCkAdd = [mCkCrd(2,:);x_crs];
                        else
                            mCkAdd = [mCkCrd(1,:);x_crs];
                        end
                    end
                    
                    if i_sgm == 1 && mTpRdi(jCrack,1)*f_xrs > s_xsg
                        
                        mTpAct(jCrack,1) = false;
                        mTpRdi(jCrack,1) = 0;
                        
                        warning(msg1,jCrack);
                        
                    end
                    
                    if i_sgm+1 == nCkCrd && mTpRdi(jCrack,2)*f_xrs > s_sgm-s_xsg
                        
                        mTpAct(jCrack,2) = false;
                        mTpRdi(jCrack,2) = 0;
                        
                        warning(msg2,jCrack);
                        
                    end
                    
                    if size(mCkAdd,1) > nSgAdd+1
                        mCkAdd = mCkAdd(end-nSgAdd:end,:);
                    end
                    
                    mCkJun(iCrack,1) = size(mCkAdd,1)-1;
                    cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                    
                    if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                            mCkJun(iCrack,2) + size(mCkAdd,1);
                    end
                    
                    mTpAct(iCrack,1) = false;
                    mTpRdi(iCrack,1) = 0;
                    
%                     break
                    
                end
            else
                
                warning(msg1,iCrack);
                mTpAct(iCrack,1) = false;
                mTpRdi(iCrack,1) = 0;
                
%                 break
                
            end
        end
    end
    
    if mTpAct(iCrack,2)
        
        %------------------------------------------------------------------
        % TIP #2
        %------------------------------------------------------------------
        
        jCkSgm = mCkJun(iCrack,:);
        nCkCrd = size(cCkCrd{iCrack},1);
        
        if jCkSgm(1)
            jCkSgm(1) = jCkSgm(1)+1;
        else
            jCkSgm(1) = 1;
        end
        
        mCkCrd_ref = cCkCrd{iCrack}(jCkSgm(1):nCkCrd,:);
        vTpCrd_ref = mCkCrd_ref(end,:);
        
        % tolerance for min. dist.
        tol_xrs = mTpRdi(iCrack,2)*f_xrs;
        tol_vtx = mTpRdi(iCrack,2)*f_vtx;
        
        % keep track of min-dist. intersections for iCrack
        mindist_r = zeros(nCrack,1); % radius to jCrack
        mindist_d = zeros(nCrack,2); % vector to jCrack
        mindist_i = zeros(nCrack,1); % segment of jCrack
        
        mindist_r(iCrack) = Inf;
        
        for jCrack = vCkAll(vCkAll ~= iCrack)
            
            jCkSgm = mCkJun(jCrack,:);
            nCkCrd = size(cCkCrd{jCrack},1);
            
            if jCkSgm(1)
                jCkSgm(1) = jCkSgm(1)+1;
            else
                jCkSgm(1) = 1;
            end
            
            if ~jCkSgm(2)
                jCkSgm(2) = nCkCrd;
            end
            
            mCkCrd = cCkCrd{jCrack}(jCkSgm(1):jCkSgm(2),:);
            nCkCrd = jCkSgm(2)-jCkSgm(1)+1;
            
            % get distance to crack
            [r,d_inc,i_sgm] = LevelSet2Poly_abs(vTpCrd_ref,mCkCrd);
            
            mindist_r(jCrack,1) = r;
            mindist_d(jCrack,:) = d_inc;
            mindist_i(jCrack,1) = i_sgm;
            
        end
        
        [r,i] = min(mindist_r);
        d_inc = mindist_d(i,:);
        i_sgm = mindist_i(i,:);
        
        jCrack = i; % min-dist. crack
        jCkSgm = mCkJun(jCrack,:);
        nCkCrd = size(cCkCrd{jCrack},1);
        
        if jCkSgm(1)
            jCkSgm(1) = jCkSgm(1)+1;
        else
            jCkSgm(1) = 1;
        end
        
        if ~jCkSgm(2)
            jCkSgm(2) = nCkCrd;
        end
        
        mCkCrd = cCkCrd{jCrack}(jCkSgm(1):jCkSgm(2),:);
        nCkCrd = jCkSgm(2)-jCkSgm(1)+1;
        
        if r < tol_xrs
            
            % tip segment direction
            d_tip = vTpCrd_ref - mCkCrd_ref(end-1,:);
            
            if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1)) ...
                    && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0
                
                % A ----x--> B
                
                x_crs = vTpCrd_ref + d_inc;
                
                x_ptA = mCkCrd(i_sgm,:);
                x_ptB = mCkCrd(i_sgm+1,:);
                
                d_sgm = x_ptB - x_ptA;
                d_xsg = x_crs - x_ptA;
                
                s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                
                % i_sgm is index into mCkCrd, then
                % j_sgm will be index into cCkCrd{jCrack}
                j_sgm = jCkSgm(1)+i_sgm-1;
                
                if s_xsg < tol_vtx % xs is on A
                    
                    if i_sgm ~= 1 || ~mCkJun(jCrack,1)
                        if i_sgm == 1 % xs is on tip
                            
                            mCkAdd = cCkCrd{jCrack};
                            mTpAct(jCrack,1) = false;
                            mTpRdi(jCrack,1) = 0;
                            
                        else
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                % mCkAdd = mCkCrd(i_sgm:end,:);
                                mCkAdd = cCkCrd{jCrack}(j_sgm:end,:);
                            else
                                % mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(j_sgm:-1:1,:);
                            end
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(1:nSgAdd+1,:);
                        end
                        
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                        cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                        
                        mTpAct(iCrack,2) = false;
                        mTpRdi(iCrack,2) = 0;
                        
%                         break
                        
                    end
                    
                elseif s_sgm-s_xsg < tol_vtx % xs on on B
                    
                    if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                        if i_sgm+1 == nCkCrd % xs is on tip
                            
                            mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                            mTpAct(jCrack,2) = false;
                            mTpRdi(jCrack,2) = 0;
                            
                        else
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                % mCkAdd = mCkCrd(i_sgm+1:end,:);
                                mCkAdd = cCkCrd{jCrack}(j_sgm+1:end,:);
                            else
                                % mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(j_sgm+1:-1:1,:);
                            end
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(1:nSgAdd+1,:);
                        end
                        
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                        cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                        
                        mTpAct(iCrack,2) = false;
                        mTpRdi(iCrack,2) = 0;
                        
%                         break
                        
                    end
                    
                else % xs is between A and B
                    
                    if nCkCrd > 2
                        if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                            % mCkAdd = [x_crs;mCkCrd(i_sgm+1:end,:)];
                            mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm+1:end,:)];
                        else
                            % mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                            mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm:-1:1,:)];
                        end
                    else
                        if s_sgm > 2*s_xsg
                            mCkAdd = [x_crs;mCkCrd(2,:)];
                        else
                            mCkAdd = [x_crs;mCkCrd(1,:)];
                        end
                    end
                    
                    if i_sgm == 1 && mTpRdi(jCrack,1)*f_xrs > s_xsg
                        
                        mTpAct(jCrack,1) = false;
                        mTpRdi(jCrack,1) = 0;
                        
                        warning(msg1,jCrack);
                        
                    end
                    
                    if i_sgm+1 == nCkCrd && mTpRdi(jCrack,2)*f_xrs > s_sgm-s_xsg
                        
                        mTpAct(jCrack,2) = false;
                        mTpRdi(jCrack,2) = 0;
                        
                        warning(msg2,jCrack);
                        
                    end
                    
                    if size(mCkAdd,1) > nSgAdd+1
                        mCkAdd = mCkAdd(1:nSgAdd+1,:);
                    end
                    
                    mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                    cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                    
                    mTpAct(iCrack,2) = false;
                    mTpRdi(iCrack,2) = 0;
                    
%                     break
                    
                end
            else
                
                warning(msg2,iCrack);
                mTpAct(iCrack,2) = false;
                mTpRdi(iCrack,2) = 0;
                
%                 break
                
            end
        end
    end
end

%--------------------------------------------------------------------------
% Check for intersections on domain boundaries
%--------------------------------------------------------------------------

global mNdCrd cBnNod;
nBound=length(cBnNod);

if nBound > 0
for i = find(mTpAct(:,1))'
    for j = 1:nBound
        
        [~,s,d] = LevelSet2Poly_bnd(cCkCrd{i}(1,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,1)*f_xrs
            if (cCkCrd{i}(1,:)-cCkCrd{i}(2,:))*d(:) > 0
                
                cCkCrd{i} = [(1+1e-4)*d+cCkCrd{i}(1,:);cCkCrd{i}];
                
                mTpRdi(i,1) = 0;
                mTpAct(i,1) = 0;
                
                if mCkJun(i,2)
                    mCkJun(i,2) = mCkJun(i,2)+1;
                end
                
                break
                
            end
        end
    end
end
for i = find(mTpAct(:,2))'
    for j = 1:nBound
        
        [~,s,d] = LevelSet2Poly_bnd(cCkCrd{i}(end,:),...
            mNdCrd(cBnNod{j}([1:end,1],1),:));
        
        if s < mTpRdi(i,2)*f_xrs
            if (cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:))*d(:) > 0
                
                cCkCrd{i} = [cCkCrd{i};cCkCrd{i}(end,:)+(1+1e-4)*d];
                
                mTpRdi(i,2) = 0;
                mTpAct(i,2) = 0;
                
                break
                
            end
        end
    end
end
end
