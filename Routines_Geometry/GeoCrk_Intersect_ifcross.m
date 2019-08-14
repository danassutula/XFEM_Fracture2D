function [mNoInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_ifcross ...
    (qCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs,mCk2Xs)
%
%==========================================================================
%
% OPTIONAL INPUT VAR.:
% mCk2Xs - cracks tips that are allowed to intersect some other crack even
%          though the latter crack has already intersected the former crack 
% qCk2Up - check these crack tips for mininum distance intersections
% mCk2Up - cracks that have just propagated.
%
%==========================================================================

% tol. for min. dist. to ck. vtx. 
f_vtx = 0.1; % as frac. of 'he'

% max n. bln. sgm. to add
nSgAdd = nSgBln; % 3; inf;

msgA = 'intersection has occurred very close to another crack tip; ';
msgB = 'a crack should not intersect another crack for the second time; ';
msgC = 'too short blending branch of slave crack jCrack';

msg1 = 'removing tip enrichment; i_crk = %d, i_tip = 1';
msg2 = 'removing tip enrichment; i_crk = %d, i_tip = 2';

% rel. to 'f_tip'
f_vtx = f_vtx/f_tip;

nCrack = length(cCkCrd);
mNoInc = zeros(nCrack,2);
vCkAll = 1:nCrack;

if isempty(mCk2Up)
    mCk2Up = false(nCrack,2);
elseif ~islogical(mCk2Up)
    if size(mCk2Up,2) == 2 && all(mCk2Up(:))
        q = mCk2Up; mCk2Up=false(nCrack,2);
        for i = q'  % [i_crack,i_tip]
            mCk2Up(i(1),i(2)) = true;
        end
    else
        error('mCk2Up dimensions are wrong')
    end
end

if ~isempty(qCk2Up) && islogical(qCk2Up)
    if size(qCk2Up,2) == 2
        [qi,qj] = find(qCk2Up');
        qCk2Up = [qj,qi];
    else
        error('qCk2Up dimensio ns are wrong')
    end
end

if nargin == 9
    mCk2Xs = false(nCrack,2);
end

for iCk2Up = qCk2Up'
    
    iCrack = iCk2Up(1);
    iCkTip = iCk2Up(2);
    
    if iCkTip == 1 && mTpAct(iCrack,1) && mCk2Up(iCrack,1)
        
        %------------------------------------------------------------------
        % TIP #1
        %------------------------------------------------------------------
        
        % get inc. sgm. coordinates
        mSgInc = cCkCrd{iCrack}(1:2,:);
        
        % get smg. direction
        d_inc = mSgInc(1,:)-mSgInc(2,:);
        u_inc = (d_inc/sqrt(d_inc(1)^2+d_inc(2)^2));
        
        % add tol. to inc. coordinate to check for intersection
        mSgInc(1,:) = mSgInc(1,:) + f_xrs*mTpRdi(iCrack,1)*u_inc;
        
        % set tol. for min. xrs. dist. to sgm. vtx.
        tol_vtx = mTpRdi(iCrack,1)*f_vtx;
        
        for jCrack = vCkAll(vCkAll ~= iCrack) % to be intersected
            
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
            
            % get crack intersection
            [x_crs,i_sgm] = GeoXPoly(mSgInc,mCkCrd,1);
            
            if ~isempty(x_crs)
                if mCk2Xs(iCrack,1)
                    pass = 1;
                elseif isempty(GeoXPoly(mCkCrd,cCkCrd{iCrack}(2:end,:),1))
                    pass = 1;
                else
                    pass = 0;
                end
            if pass
         
                % A ----x--> B

                i_sgm = i_sgm(1,2);

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
                        if i_sgm == 1 % xs is on crack tip
                            
                            if mCk2Up(jCrack,2)
                                % incremented OR intersecting
                                mCkAdd = mCkCrd(end-1:-1:1,:);
                            else
                                % mCkAdd = mCkCrd(end:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                            end

                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;

                        else % i_sgm ~= 1
                            
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(end-1:-1:i_sgm,:); 
                                else
                                    % mCkAdd = mCkCrd(end:-1:i_sgm,:);
                                    mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,1)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(2:i_sgm,:);
                                    else
                                        % mCkAdd = mCkCrd(1:i_sgm,:);
                                        mCkAdd = cCkCrd{jCrack}(1:j_sgm,:);
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(2:i_sgm,:);
                                else
                                    % mCkAdd = mCkCrd(1:i_sgm,:);
                                    mCkAdd = cCkCrd{jCrack}(1:j_sgm,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,2)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(end-1:-1:i_sgm,:);
                                    else
                                        % mCkAdd = mCkCrd(end:-1:i_sgm,:);
                                        mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm,:);
                                    end
                                end
                                
                            end
                        end
                        
                        if size(mCkAdd,1) < 2
                            warning(msgC)
                            keyboard
                        end
                        
                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(end-nSgAdd:end,:);
                        end

                        mNoInc(iCrack,1) = size(mCkAdd,1)-1;
                        mCkJun(iCrack,1) = mNoInc(iCrack,1);
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];

                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                            mCkJun(iCrack,2) + mNoInc(iCrack,1);
                        end

                        mTpAct(iCrack,1) = false;
                        mTpRdi(iCrack,1) = 0;

                        break

                    end

                elseif s_sgm-s_xsg < tol_vtx % xs on B

                    if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                        if i_sgm+1 == nCkCrd % xs is on crack tip
                            
                            if mCk2Up(jCrack,1)
                                % incremented OR intersecting
                                mCkAdd = mCkCrd(2:end,:);
                            else
                                % mCkAdd = mCkCrd;
                                mCkAdd = cCkCrd{jCrack};
                            end
                            
                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;

                        else % i_sgm+1 ~= nCkCrd
                            
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(end-1:-1:i_sgm+1,:);
                                else
                                    % mCkAdd = mCkCrd(end:-1:i_sgm+1,:);
                                    mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm+1,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,1)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(2:i_sgm+1,:);
                                    else
                                        % mCkAdd = mCkCrd(1:i_sgm+1,:);
                                        mCkAdd = cCkCrd{jCrack}(1:j_sgm+1,:);
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(2:i_sgm+1,:);
                                else
                                    % mCkAdd = mCkCrd(1:i_sgm+1,:);
                                    mCkAdd = cCkCrd{jCrack}(1:j_sgm+1,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,2)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(end-1:-1:i_sgm+1,:);
                                    else
                                        % mCkAdd = mCkCrd(end:-1:i_sgm+1,:);
                                        mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm+1,:);
                                    end
                                end
                                
                            end
                        end

                        if size(mCkAdd,1) < 2
                            warning(msgC)
                            keyboard
                        end

                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(end-nSgAdd:end,:);
                        end

                        mNoInc(iCrack,1) = size(mCkAdd,1)-1;
                        mCkJun(iCrack,1) = mNoInc(iCrack,1);
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];

                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                            mCkJun(iCrack,2) + mNoInc(iCrack,1);
                        end
                        
                        mTpAct(iCrack,1) = false;
                        mTpRdi(iCrack,1) = 0;

                        break

                    end

                else % xs is between A and B
                    
                    if nCkCrd > 2
                        if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                            
                            if mCk2Up(jCrack,2)
                                % incremented OR intersecting
                                mCkAdd = [mCkCrd(end-1:-1:i_sgm+1,:);x_crs];
                            else
                                % mCkAdd = [mCkCrd(end:-1:i_sgm+1,:);x_crs];
                                mCkAdd = [cCkCrd{jCrack}(end:-1:j_sgm+1,:);x_crs];
                            end
                            if size(mCkAdd,1)<2
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = [mCkCrd(2:i_sgm,:);x_crs];
                                else
                                    % mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                                    mCkAdd = [cCkCrd{jCrack}(1:j_sgm,:);x_crs];
                                end
                            end
                            
                        else
                            
                            if mCk2Up(jCrack,1)
                                % incremented OR intersecting
                                mCkAdd = [mCkCrd(2:i_sgm,:);x_crs];
                            else
                                % mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                                mCkAdd = [cCkCrd{jCrack}(1:j_sgm,:);x_crs];
                            end
                            if size(mCkAdd,1)<2
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = [mCkCrd(end-1:-1:i_sgm+1,:);x_crs];
                                else
                                    % mCkAdd = [mCkCrd(end:-1:i_sgm+1,:);x_crs];
                                    mCkAdd = [cCkCrd{jCrack}(end:-1:j_sgm+1,:);x_crs];
                                end
                            end
                            
                        end
                    else
                        if s_sgm > 2*s_xsg
                            mCkAdd = [mCkCrd(2,:);x_crs];
                        else
                            mCkAdd = [mCkCrd(1,:);x_crs];
                        end
                    end

                    if size(mCkAdd,1) < 2
                        warning(msgC)
                        keyboard
                    end

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
                    end
                    
                    mNoInc(iCrack,1) = size(mCkAdd,1)-1;
                    mCkJun(iCrack,1) = mNoInc(iCrack,1);
                    cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}(2:end,:)];
                    
                    if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                        mCkJun(iCrack,2) + mNoInc(iCrack,1);
                    end
                    
                    mTpAct(iCrack,1) = false;
                    mTpRdi(iCrack,1) = 0;

                    break

                end

            else

                cCkCrd{iCrack} = cCkCrd{iCrack}(2:end,:);

                mTpAct(iCrack,1) = false;
                mTpRdi(iCrack,1) = 0;
                mNoInc(iCrack,1) =-1;

                warning([msgB,msg1],iCrack);
                break

            end
            end
        end
        
    elseif iCkTip == 2 && mTpAct(iCrack,2) && mCk2Up(iCrack,2)
        
        %------------------------------------------------------------------
        % TIP #2
        %------------------------------------------------------------------
        
        % get inc. sgm. coordinates
        mSgInc = cCkCrd{iCrack}(end-1:end,:);
        
        % get inc. direction
        d_inc = mSgInc(2,:)-mSgInc(1,:);
        u_inc = d_inc/sqrt(d_inc(1)^2+d_inc(2)^2);
        
        % add tol. to inc. coordinate to check for intersection
        mSgInc(2,:) = mSgInc(2,:) + f_xrs*mTpRdi(iCrack,2)*u_inc;
        
        % set tol. for min. xrs. dist. to sgm. vtx.
        tol_vtx = mTpRdi(iCrack,2)*f_vtx;
        
        for jCrack = vCkAll(vCkAll ~= iCrack) % to be intersected
            
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
            
            % get crack intersection
            [x_crs,i_sgm] = GeoXPoly(mSgInc,mCkCrd,1);
            
            if ~isempty(x_crs)
                if mCk2Xs(iCrack,2)
                    pass = 1;
                elseif isempty(GeoXPoly(mCkCrd,cCkCrd{iCrack}(1:end-1,:),1))
                    pass = 1;
                else
                    pass = 0;
                end
            if pass
                
                % A ----x--> B

                i_sgm = i_sgm(1,2);

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
                        if i_sgm == 1 % xs is on crack tip
                            
                            if mCk2Up(jCrack,2)
                                % incremented OR intersecting
                                mCkAdd = mCkCrd(1:end-1,:);
                            else
                                % mCkAdd = mCkCrd;
                                mCkAdd = cCkCrd{jCrack};
                            end

                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;

                        else % i_sgm ~= 1
                            
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(i_sgm:end-1,:);
                                else
                                    % mCkAdd = mCkCrd(i_sgm:end,:);
                                    mCkAdd = cCkCrd{jCrack}(j_sgm:end,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,1)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(i_sgm:-1:2,:);
                                    else
                                        % mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm:-1:1,:);
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(i_sgm:-1:2,:);
                                else
                                    % mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                    mCkAdd = cCkCrd{jCrack}(j_sgm:-1:1,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,2)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(i_sgm:end-1,:);
                                    else
                                        % mCkAdd = mCkCrd(i_sgm:end,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm:end,:);
                                    end
                                end
                                
                            end
                        end

                        if size(mCkAdd,1) < 2
                            warning(msgC)
                            keyboard
                        end

                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(1:nSgAdd+1,:);
                        end
                        
                        mNoInc(iCrack,2) = size(mCkAdd,1)-1;
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack});
                        cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];
                        
                        mTpAct(iCrack,2) = false;
                        mTpRdi(iCrack,2) = 0;

                        break

                    end

                elseif s_sgm-s_xsg < tol_vtx % xs on on B

                    if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                        if i_sgm+1 == nCkCrd % xs is on crack tip
                            
                            if mCk2Up(jCrack,1)
                                % incremented OR intersecting
                                mCkAdd = mCkCrd(end:-1:2,:);
                            else
                                % mCkAdd = mCkCrd(end:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                            end

                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;

                        else % i_sgm+1 ~= nCkCrd
                            
                            if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(i_sgm+1:end-1,:);
                                else
                                    % mCkAdd = mCkCrd(i_sgm+1:end,:);
                                    mCkAdd = cCkCrd{jCrack}(j_sgm+1:end,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,1)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(i_sgm+1:-1:2,:);
                                    else
                                        % mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm+1:-1:1,:);
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = mCkCrd(i_sgm+1:-1:2,:);
                                else
                                    % mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                    mCkAdd = cCkCrd{jCrack}(j_sgm+1:-1:1,:);
                                end
                                if size(mCkAdd,1)<2
                                    if mCk2Up(jCrack,2)
                                        % incremented OR intersecting
                                        mCkAdd = mCkCrd(i_sgm+1:end-1,:);
                                    else
                                        % mCkAdd = mCkCrd(i_sgm+1:end,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm+1:end,:);
                                    end
                                end
                                
                            end
                        end
                        
                        if size(mCkAdd,1) < 2
                            warning(msgC)
                            keyboard
                        end

                        if size(mCkAdd,1) > nSgAdd+1
                            mCkAdd = mCkAdd(1:nSgAdd+1,:);
                        end
                        
                        mNoInc(iCrack,2) = size(mCkAdd,1)-1;
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack});
                        cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];

                        mTpAct(iCrack,2) = false;
                        mTpRdi(iCrack,2) = 0;

                        break

                    end

                else % xs is between A and B

                    if nCkCrd > 2
                        if d_inc(1)*d_sgm(1)+d_inc(2)*d_sgm(2) > 0
                            
                            if mCk2Up(jCrack,2)
                                % incremented OR intersecting
                                mCkAdd = [x_crs;mCkCrd(i_sgm+1:end-1,:)];
                            else
                                % mCkAdd = [x_crs;mCkCrd(i_sgm+1:end,:)];
                                mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm+1:end,:)];
                            end
                            if size(mCkAdd,1)<2
                                if mCk2Up(jCrack,1)
                                    % incremented OR intersecting
                                    mCkAdd = [x_crs;mCkCrd(i_sgm:-1:2,:)];
                                else
                                    % mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                                    mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm:-1:1,:)];
                                end
                            end
                            
                        else
                            
                            if mCk2Up(jCrack,1)
                                % incremented OR intersecting
                                mCkAdd = [x_crs;mCkCrd(i_sgm:-1:2,:)];
                            else
                                % mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                                mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm:-1:1,:)];
                            end
                            if size(mCkAdd,1)<2
                                if mCk2Up(jCrack,2)
                                    % incremented OR intersecting
                                    mCkAdd = [x_crs;mCkCrd(i_sgm+1:end-1,:)];
                                else
                                    % mCkAdd = [x_crs;mCkCrd(i_sgm+1:end,:)];
                                    mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm+1:end,:)];
                                end
                            end
                            
                        end
                    else
                        if s_sgm > 2*s_xsg
                            mCkAdd = [x_crs;mCkCrd(2,:)];
                        else
                            mCkAdd = [x_crs;mCkCrd(1,:)];
                        end
                    end
                    
                    if size(mCkAdd,1) < 2
                        warning(msgC)
                        keyboard
                    end
                    
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
                        mCkAdd = mCkAdd(1:nSgAdd+1,:);
                    end
                    
                    mNoInc(iCrack,2) = size(mCkAdd,1)-1;
                    mCkJun(iCrack,2) = length(cCkCrd{iCrack});
                    cCkCrd{iCrack} = [cCkCrd{iCrack}(1:end-1,:);mCkAdd];

                    mTpAct(iCrack,2) = false;
                    mTpRdi(iCrack,2) = 0;

                    break

                end
                
            else

                cCkCrd{iCrack} = cCkCrd{iCrack}(1:end-1,:);

                mTpAct(iCrack,2) = false;
                mTpRdi(iCrack,2) = 0;
                mNoInc(iCrack,2) =-1;

                warning([msgB,msg2],iCrack);
                break

            end
            end
        end
    end
end