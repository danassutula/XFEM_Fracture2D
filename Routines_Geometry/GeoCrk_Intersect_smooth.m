function [mNoInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_smooth...
    (qCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs)

error('Where is this script?')

%==========================================================================

% mCk2Ch - check these crack tips for min. dist. tip-to-tip intersection
% mCk2Up - cracks that have just propagated.

%==========================================================================

% tol. for min. dist. to ck. vtx. 
f_vtx = 0.1; % as frac. of 'he'

% max n. bln. sgm. to add
nSgAdd = nSgBln; % 3; inf;

msgA = 'intersection has occurred very close to another crack tip; ';
msgB = 'a crack should not intersect another crack for the second time; ';

msg1 = 'removing tip enrichment; i_crk = %d, i_tip = 1';
msg2 = 'removing tip enrichment; i_crk = %d, i_tip = 2';

% rel. to 'f_tip'
f_vtx = f_vtx/f_tip;

nCrack = length(cCkCrd);
mNoInc = zeros(nCrack,2);
vCkAll = 1:nCrack;

if isempty(mCk2Up)
    mCk2Up = true(nCrack,2);
end

if isempty(qCk2Up)
    [qi,qj] = find(mCk2Up');
    qCk2Up = [qj,qi];
end

% iCrack - the intersecting crack
% jCrack - the intersected crack

for iCk2Up = qCk2Up'
    
    iCrack = iCk2Up(1);
    iCkTip = iCk2Up(2);
    
    if iCkTip == 1 && mTpAct(iCrack,1) && mCk2Up(iCrack,1)
        
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
            
            if r < tol_xrs
                
                % tip segment direction
                d_tip = vTpCrd_ref - mCkCrd_ref(2,:);
                
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1)) % ...
                    % && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0
                    
                    % A ----x--> B
                    
                    x_crs = vTpCrd_ref + d_inc;
                    
                    x_ptA = mCkCrd(i_sgm,:);
                    x_ptB = mCkCrd(i_sgm+1,:);
                    
                    d_sgm = x_ptB - x_ptA;
                    d_xsg = x_crs - x_ptA;
                    
                    s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                    s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                    
                    % modify inc. direction
                    d_inc = d_inc + d_tip;
                    
                    if s_xsg < tol_vtx % xs is on A
                        
                        if i_sgm ~= 1 || ~mCkJun(jCrack,1)
                            
                            if i_sgm == 1 % xs is on crack tip
                                
                                if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                    mCkAdd = mCkCrd(nCkCrd-1:-1:1,:);
                                    mCkJun(iCrack,1) = nCkCrd-2;
                                else
                                    mCkAdd = mCkCrd(nCkCrd:-1:1,:);
                                    mCkJun(iCrack,1) = nCkCrd-1;
                                end
                                
                                mTpAct(jCrack,1) = false;
                                mCk2Up(jCrack,1) = true;
                                mTpRdi(jCrack,1) = 0;
                                
                            else
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                        mCkAdd = mCkCrd(nCkCrd-1:-1:i_sgm,:);
                                        mCkJun(iCrack,1) = nCkCrd-1-i_sgm;
                                    else
                                        mCkAdd = mCkCrd(nCkCrd:-1:i_sgm,:);
                                        mCkJun(iCrack,1) = nCkCrd-i_sgm;
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                        mCkAdd = mCkCrd(2:i_sgm,:);
                                        mCkJun(iCrack,1) = i_sgm-2;
                                    else
                                        mCkAdd = mCkCrd(1:i_sgm,:);
                                        mCkJun(iCrack,1) = i_sgm-1;
                                    end
                                    
                                end
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(end-nSgAdd:end,:);
                                mCkJun(iCrack,1) = nSgAdd;
                            end
                            
                            if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                                mCkJun(iCrack,2) + size(mCkAdd,1)-1;
                            end
                            
                            cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];
                            
                            mNoInc(iCrack,1) = size(mCkAdd,1) - 1;
                            mTpAct(iCrack,1) = false;
                            mTpRdi(iCrack,1) = 0;
                            
                            break
                            
                        end
                        
                    elseif s_sgm-s_xsg < tol_vtx % xs on B
                        
                        if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                            
                            if i_sgm+1 == nCkCrd % xs is on crack tip
                                
                                if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                    mCkAdd = mCkCrd(2:nCkCrd,:);
                                    mCkJun(iCrack,1) = i_sgm-1;
                                else
                                    mCkAdd = mCkCrd(1:nCkCrd,:);
                                    mCkJun(iCrack,1) = i_sgm;
                                end
                                
                                mTpAct(jCrack,2) = false;
                                mCk2Up(jCrack,2) = true;
                                mTpRdi(jCrack,2) = 0;
                                
                            else
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                        mCkAdd = mCkCrd(nCkCrd-1:-1:i_sgm+1,:);
                                        mCkJun(iCrack,1) = nCkCrd-1-i_sgm-1;
                                    else
                                        mCkAdd = mCkCrd(nCkCrd:-1:i_sgm+1,:);
                                        mCkJun(iCrack,1) = nCkCrd-i_sgm-1;
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                        mCkAdd = mCkCrd(2:i_sgm+1,:);
                                        mCkJun(iCrack,1) = i_sgm-1;
                                    else 
                                        mCkAdd = mCkCrd(1:i_sgm+1,:);
                                        mCkJun(iCrack,1) = i_sgm;
                                    end
                                    
                                end
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(end-nSgAdd:end,:);
                                mCkJun(iCrack,1) = nSgAdd;
                            end
                            
                            if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                                mCkJun(iCrack,2) + size(mCkAdd,1)-1;
                            end
                            
                            cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];
                            
                            mNoInc(iCrack,1) = size(mCkAdd,1) - 1;
                            mTpAct(iCrack,1) = false;
                            mTpRdi(iCrack,1) = 0;
                            
                            break
                            
                        end
                        
                    else % xs is between A and B
                        
%                         if nCkCrd > 2
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                    mCkAdd = [mCkCrd(nCkCrd-1:-1:i_sgm+1,:);x_crs];
                                    mCkJun(iCrack,1) = nCkCrd-1-i_sgm;
                                else
                                    mCkAdd = [mCkCrd(nCkCrd:-1:i_sgm+1,:);x_crs];
                                    mCkJun(iCrack,1) = nCkCrd-i_sgm;
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                    mCkAdd = [mCkCrd(2:i_sgm,:);x_crs];
                                    mCkJun(iCrack,1) = i_sgm-1;
                                else
                                    mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                                    mCkJun(iCrack,1) = i_sgm;
                                end
                                
                            end
%                         else
%                             if s_sgm > 2*s_xsg
%                                 mCkAdd = [mCkCrd(2,:);x_crs];
%                                 mCkJun(iCrack,1) = 1;
%                             else
%                                 mCkAdd = [mCkCrd(1,:);x_crs];
%                                 mCkJun(iCrack,1) = 1;
%                             end
%                         end
                        
                        if i_sgm == 1 && mTpRdi(jCrack,1)*f_xrs > s_xsg
                            
                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;
                            
                            warning([msgA,msg1],jCrack);
                            
                        end
                        if i_sgm+1 == nCkCrd && mTpRdi(jCrack,2)*f_xrs > s_sgm-s_xsg
                            
                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;
                            
                            warning([msgA,msg2],jCrack);
                            
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(end-nSgAdd:end,:);
                            mCkJun(iCrack,1) = nSgAdd;
                        end
                        
                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                            mCkJun(iCrack,2) + size(mCkAdd,1)-1;
                        end
                        
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];
                        
                        mNoInc(iCrack,1) = size(mCkAdd,1) - 1;
                        mTpAct(iCrack,1) = false;
                        mTpRdi(iCrack,1) = 0;
                        
                        break
                        
                    end
                else
                    
                    warning([msgB,msg1],iCrack);
                    mTpAct(iCrack,1) = false;
                    mTpRdi(iCrack,1) = 0;
                    
                    break
                    
                end
            end
        end
    
    elseif iCkTip == 2 && mTpAct(iCrack,2) && mCk2Up(iCrack,2)
        
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
            
            if r < tol_xrs
                
                % tip segment direction
                d_tip = vTpCrd_ref - mCkCrd_ref(end-1,:);
                
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1)) % ...
                    % && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0
                    
                    % A ----x--> B
                    
                    x_crs = vTpCrd_ref + d_inc;
                    
                    x_ptA = mCkCrd(i_sgm,:);
                    x_ptB = mCkCrd(i_sgm+1,:);
                    
                    d_sgm = x_ptB - x_ptA;
                    d_xsg = x_crs - x_ptA;
                    
                    s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                    s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                    
                    % modify inc. direction
                    d_inc = d_inc + d_tip;
                    
                    if s_xsg < tol_vtx % xs is on A
                        
                        if i_sgm ~= 1 || ~mCkJun(jCrack,1)
                            
                            if i_sgm == 1 % xs is on crack tip
                                
                                if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                    mCkAdd = mCkCrd(1:nCkCrd-1,:);
                                else
                                    mCkAdd = mCkCrd(1:nCkCrd,:);
                                end
                                
                                mTpAct(jCrack,1) = false;
                                mCk2Up(jCrack,1) = true;
                                mTpRdi(jCrack,1) = 0;
                                
                            else
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                        mCkAdd = mCkCrd(i_sgm:nCkCrd-1,:);
                                    else
                                        mCkAdd = mCkCrd(i_sgm:nCkCrd,:);
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                        mCkAdd = mCkCrd(i_sgm:-1:2,:);
                                    else
                                        mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                    end
                                    
                                end
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(1:nSgAdd+1,:);
                            end
                            
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack}); % -1 + 1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1) - 1;
                            mTpAct(iCrack,2) = false;
                            mTpRdi(iCrack,2) = 0;
                            
                            break
                            
                        end
                        
                    elseif s_sgm-s_xsg < tol_vtx % xs on on B
                        
                        if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                            
                            if i_sgm+1 == nCkCrd % xs is on crack tip
                                
                                if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                    mCkAdd = mCkCrd(nCkCrd:-1:2,:);
                                else
                                    mCkAdd = mCkCrd(nCkCrd:-1:1,:);
                                end
                                
                                mTpAct(jCrack,2) = false;
                                mCk2Up(jCrack,2) = true;
                                mTpRdi(jCrack,2) = 0;
                                
                            else
                                
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                        mCkAdd = mCkCrd(i_sgm+1:nCkCrd-1,:);
                                    else
                                        mCkAdd = mCkCrd(i_sgm+1:nCkCrd,:);
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                        mCkAdd = mCkCrd(i_sgm+1:-1:2,:);
                                    else
                                        mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                    end
                                    
                                end
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(1:nSgAdd+1,:);
                            end
                            
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack}); % -1 + 1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1) - 1;
                            mTpAct(iCrack,2) = false;
                            mTpRdi(iCrack,2) = 0;
                            
                            break
                            
                        end
                        
                    else % xs is between A and B
                        
%                         if nCkCrd > 2
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2) && mTpAct(jCrack,2) % incremented tip
                                    mCkAdd = [x_crs;mCkCrd(i_sgm+1:end-1,:)];
                                else
                                    mCkAdd = [x_crs;mCkCrd(i_sgm+1:end,:)];
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1) && mTpAct(jCrack,1) % incremented tip
                                    mCkAdd = [x_crs;mCkCrd(i_sgm:-1:2,:)];
                                else
                                    mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                                end
                                
                            end
%                         else
%                             if s_sgm > 2*s_xsg
%                                 mCkAdd = [x_crs;mCkCrd(2,:)];
%                             else
%                                 mCkAdd = [x_crs;mCkCrd(1,:)];
%                             end
%                         end
                        
                        if i_sgm == 1 && mTpRdi(jCrack,1)*f_xrs > s_xsg
                            
                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;
                            
                            warning([msgA,msg1],jCrack);
                            
                        elseif i_sgm+1 == nCkCrd && mTpRdi(jCrack,2)*f_xrs > s_sgm-s_xsg
                            
                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;
                            
                            warning([msgA,msg2],jCrack);
                            
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(1:nSgAdd+1,:);
                        end
                        
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack}); % -1 + 1;
                        cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];
                        
                        mNoInc(iCrack,2) = size(mCkAdd,1) - 1;
                        mTpAct(iCrack,2) = false;
                        mTpRdi(iCrack,2) = 0;
                        
                        break
                        
                    end
                else
                    
                    warning([msgB,msg2],iCrack);
                    mTpAct(iCrack,2) = false;
                    mTpRdi(iCrack,2) = 0;
                    
                    break
                    
                end
            end
        end
    end
end