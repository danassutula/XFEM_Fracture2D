function [mNoInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_mindist ...
    (qCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs)

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
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1))
                    
                    % tip segment direction
                    d_tip = vTpCrd_ref - mCkCrd_ref(2,:);
                    
                    % if kink angle is obtuse, make it pi/2
                    if d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) < 0
                        
                        d_nml = [-d_tip(2),d_tip(1)];
                        d_nml = d_nml*sign(d_nml*d_inc(:));
                        
                        [x_crs,i_sgx] = GeoXPoly([vTpCrd_ref;...
                            vTpCrd_ref + f_xrs*d_nml],mCkCrd,1);
                        
                        if size(x_crs,1) == 1
                            % make inc. normal to tip
                            d_inc = x_crs-vTpCrd_ref;
                            i_sgm = i_sgx(2);
                        end
                        
                    end
                    
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
                            if i_sgm == 1 % xs is on crack tip
                                
                                if mCk2Up(jCrack,2) % incremented tip
                                    if mCkJun(jCrack,2) && nCkCrd>3
                                        mCkAdd = mCkCrd(end-2:-1:1,:);
                                        % must be that: size(mCkCrd,1)>=2
                                    else
                                        mCkAdd = mCkCrd(end-1:-1:1,:);
                                        % must be that: size(mCkCrd,1)>=2
                                    end
                                else
                                    % mCkAdd = mCkCrd(end:-1:1,:);
                                    mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                                    % must be that: size(mCkCrd,1)>=2
                                end
                                
                                mTpAct(jCrack,1) = false;
                                mCk2Up(jCrack,1) = true;
                                mTpRdi(jCrack,1) = 0;
                                
                                % will always be sufficient:
                                % must be that size(mCkCrd,1)>=2
                                
                            else
                                
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCrack-i_sgm>2
                                            mCkAdd = mCkCrd(end-2:-1:i_sgm,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = mCkCrd(end-1:-1:i_sgm,:);
                                            % possibly: size(mCkCrd,1)<2
                                        end
                                    else
                                        % mCkAdd = mCkCrd(end:-1:i_sgm,:);
                                        mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm,:);
                                        % possibly: size(mCkCrd,1)<2
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,1)
                                            if mCkJun(jCrack,1) && i_sgm>3
                                                mCkAdd = mCkCrd(3:i_sgm,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            else
                                                mCkAdd = mCkCrd(2:i_sgm,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            end
                                        else
                                            mCkAdd = cCkCrd{jCrack}(1:j_sgm,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>3
                                            mCkAdd = mCkCrd(3:i_sgm,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = mCkCrd(2:i_sgm,:);
                                            % possibly: size(mCkCrd,1)<2
                                        end
                                    else
                                        % mCkAdd = mCkCrd(1:i_sgm,:);
                                        mCkAdd = cCkCrd{jCrack}(1:j_sgm,:);
                                        % possibly: size(mCkCrd,1)<2
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,2)
                                            if mCkJun(jCrack,2) && nCrack-i_sgm>2
                                                mCkAdd = mCkCrd(end-2:-1:i_sgm,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            else
                                                mCkAdd = mCkCrd(end-1:-1:i_sgm,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            end
                                        else
                                            mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    end
                                end
                            end
                             
                            if size(mCkAdd,1) < 2
                                warning(msgC)
                                keyboard
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                % usually, unnecessary blending length
                                mCkAdd = mCkAdd(end-nSgAdd:end,:);
                            end
                            
                            mNoInc(iCrack,1) = size(mCkAdd,1);
                            mCkJun(iCrack,1) = mNoInc(iCrack,1)-1;
                            cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                            
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
                                
                                if mCk2Up(jCrack,1) % just incremented
                                    if mCkJun(jCrack,1) && nCkCrd>3
                                        mCkAdd = mCkCrd(3:end,:);
                                        % must be that: size(mCkCrd,1)>=2
                                    else
                                        mCkAdd = mCkCrd(2:end,:);
                                    end
                                else
                                    % mCkAdd = mCkCrd;
                                    mCkAdd = cCkCrd{jCrack};
                                end
                                
                                mTpAct(jCrack,2) = false;
                                mCk2Up(jCrack,2) = true;
                                mTpRdi(jCrack,2) = 0;
                                
                                % will always be sufficient:
                                % must be that size(mCkCrd,1)>=2
                                
                            else
                                
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCkCrd-i_sgm>3
                                            mCkAdd = mCkCrd(end-2:-1:i_sgm+1,:);
                                            % must be that size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = mCkCrd(end-1:-1:i_sgm+1,:);
                                            % possibly size(mCkCrd,1)<2
                                        end
                                    else
                                        % mCkAdd = mCkCrd(end:-1:i_sgm+1,:);
                                        mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm+1,:);
                                        % possibly size(mCkCrd,1)<2
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,1) % just incremented
                                            if mCkJun(jCrack,1) && i_sgm>2
                                                mCkAdd = mCkCrd(3:i_sgm+1,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            else
                                                mCkAdd = mCkCrd(2:i_sgm+1,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            end
                                        else
                                            mCkAdd = cCkCrd{jCrack}(1:j_sgm+1,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>2
                                            mCkAdd = mCkCrd(3:i_sgm+1,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = mCkCrd(2:i_sgm+1,:);
                                            % possibly size(mCkCrd,1)<2
                                        end
                                    else 
                                        % mCkAdd = mCkCrd(1:i_sgm+1,:);
                                        mCkAdd = cCkCrd{jCrack}(1:j_sgm+1,:);
                                        % possibly size(mCkCrd,1)<2
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,2) % just incremented
                                            if mCkJun(jCrack,2) && nCkCrd-i_sgm>3
                                                mCkAdd = mCkCrd(end-2:-1:i_sgm+1,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            else
                                                mCkAdd = mCkCrd(end-1:-1:i_sgm+1,:);
                                                % must be that: size(mCkCrd,1)>=2
                                            end
                                        else
                                            % mCkAdd = mCkCrd(end:-1:i_sgm+1,:);
                                            mCkAdd = cCkCrd{jCrack}(end:-1:j_sgm+1,:);
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    end
                                end
                            end
                            
                            if size(mCkAdd,1) < 2
                                warning(msgC)
                                keyboard
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                % usually, unnecessary blending length
                                mCkAdd = mCkAdd(end-nSgAdd:end,:);
                            end
                            
                            mNoInc(iCrack,1) = size(mCkAdd,1);
                            mCkJun(iCrack,1) = mNoInc(iCrack,1)-1;
                            cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                            
                            if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                                mCkJun(iCrack,2) + mNoInc(iCrack,1);
                            end
                            
                            mTpAct(iCrack,1) = false;
                            mTpRdi(iCrack,1) = 0;
                            
                            break
                            
                        end
                        
                    else % xs is between A and B
                        
                        if nCkCrd > 2
                            if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2) % just incremented
                                    if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                        mCkAdd = [mCkCrd(end-2:-1:i_sgm+1,:);x_crs];
                                        % must be that: size(mCkCrd,1)>=2
                                    else
                                        mCkAdd = [mCkCrd(end-1:-1:i_sgm+1,:);x_crs];
                                        % possibly: size(mCkCrd,1)<2
                                    end
                                else
                                    % mCkAdd = [mCkCrd(end:-1:i_sgm+1,:);x_crs];
                                    mCkAdd = [cCkCrd{jCrack}(end:-1:j_sgm+1,:);x_crs];
                                    % possibly: size(mCkCrd,1)<2
                                end
                                if size(mCkAdd,1)<2
                                    % inadequate slave crack deflection direction
                                    % try deflecting jCrack in opposite direction
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>2
                                            mCkAdd = [mCkCrd(3:i_sgm,:);x_crs];
                                            % must be that: size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = [mCkCrd(2:i_sgm,:);x_crs];
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    else
                                        % mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                                        mCkAdd = [cCkCrd{jCrack}(1:j_sgm,:);x_crs];
                                        % must be that: size(mCkCrd,1)>=2
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1) % just incremented
                                    if mCkJun(jCrack,1) && i_sgm>2
                                        mCkAdd = [mCkCrd(3:i_sgm,:);x_crs];
                                        % must be that: size(mCkCrd,1)>=2
                                    else
                                        mCkAdd = [mCkCrd(2:i_sgm,:);x_crs];
                                        % possibly: size(mCkCrd,1)<2
                                    end
                                else
                                    % mCkAdd = [mCkCrd(1:i_sgm,:);x_crs];
                                    mCkAdd = [cCkCrd{jCrack}(1:j_sgm,:);x_crs];
                                end
                                if size(mCkAdd,1)<2
                                    % inadequate slave crack deflection direction
                                    % try deflecting jCrack in opposite direction
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                            mCkAdd = [mCkCrd(end-2:-1:i_sgm+1,:);x_crs];
                                            % must be that: size(mCkCrd,1)>=2
                                        else
                                            mCkAdd = [mCkCrd(end-1:-1:i_sgm+1,:);x_crs];
                                            % must be that: size(mCkCrd,1)>=2
                                        end
                                    else
                                        % mCkAdd = [mCkCrd(end:-1:i_sgm+1,:);x_crs];
                                        mCkAdd = [cCkCrd{jCrack}(end:-1:j_sgm+1,:);x_crs];
                                        % must be that: size(mCkCrd,1)>=2
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
                        
                        mNoInc(iCrack,1) = size(mCkAdd,1);
                        mCkJun(iCrack,1) = mNoInc(iCrack,1)-1;
                        cCkCrd{iCrack} = [mCkAdd;cCkCrd{iCrack}];
                        
                        if mCkJun(iCrack,2); mCkJun(iCrack,2) = ...
                            mCkJun(iCrack,2) + mNoInc(iCrack,1);
                        end
                        
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
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1))
                    
                    % tip segment direction
                    d_tip = vTpCrd_ref - mCkCrd_ref(end-1,:);
                    
                    % if kink angle is obtuse, make it pi/2
                    if d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) < 0
                        
                        d_nml = [-d_tip(2),d_tip(1)];
                        d_nml = d_nml*sign(d_nml*d_inc(:));
                        
                        [x_crs,i_sgx] = GeoXPoly([vTpCrd_ref;...
                            vTpCrd_ref + f_xrs*d_nml],mCkCrd,1);
                        
                        if size(x_crs,1) == 1
                            % make inc. normal to tip
                            d_inc = x_crs-vTpCrd_ref;
                            i_sgm = i_sgx(2);
                        end
                        
                    end
                    
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
                                
                                if mCk2Up(jCrack,2) % just incremented
                                    if mCkJun(jCrack,2) && nCkCrd>3
                                        mCkAdd = mCkCrd(1:end-2,:);
                                    else
                                        mCkAdd = mCkCrd(1:end-1,:);
                                    end
                                else    
                                    % mCkAdd = mCkCrd;
                                    mCkAdd = cCkCrd{jCrack};
                                end
                                
                                mTpAct(jCrack,1) = false;
                                mCk2Up(jCrack,1) = true;
                                mTpRdi(jCrack,1) = 0;
                                
                                % will always be sufficient:
                                % must be that size(mCkCrd,1)>=2
                                
                            else
                                
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                            mCkAdd = mCkCrd(i_sgm:end-2,:);
                                        else
                                            mCkAdd = mCkCrd(i_sgm:end-1,:);
                                        end
                                    else
                                        % mCkAdd = mCkCrd(i_sgm:end,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm:end,:);
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,1) % just incremented
                                            if mCkJun(jCrack,1) && i_sgm>3
                                                mCkAdd = mCkCrd(i_sgm:-1:3,:);
                                            else
                                                mCkAdd = mCkCrd(i_sgm:-1:2,:);
                                            end
                                        else
                                            % mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                            mCkAdd = cCkCrd{jCrack}(j_sgm:-1:1,:);
                                        end
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>3
                                            mCkAdd = mCkCrd(i_sgm:-1:3,:);
                                        else
                                            mCkAdd = mCkCrd(i_sgm:-1:2,:);
                                        end
                                    else
                                        % mCkAdd = mCkCrd(i_sgm:-1:1,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm:-1:1,:);
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,2) % just incremented
                                            if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                                mCkAdd = mCkCrd(i_sgm:end-2,:);
                                            else
                                                mCkAdd = mCkCrd(i_sgm:end-1,:);
                                            end
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
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1);
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                            
                            mTpAct(iCrack,2) = false;
                            mTpRdi(iCrack,2) = 0;
                            
                            break
                            
                        end
                        
                    elseif s_sgm-s_xsg < tol_vtx % xs on on B
                        
                        if i_sgm+1 ~= nCkCrd || ~mCkJun(jCrack,2)
                            if i_sgm+1 == nCkCrd % xs is on tip
                                
                                if mCk2Up(jCrack,1) % just incremented
                                    if mCkJun(jCrack,1) && nCkCrd>3
                                        mCkAdd = mCkCrd(end:-1:3,:);
                                    else
                                        mCkAdd = mCkCrd(end:-1:2,:);
                                    end
                                else
                                    % mCkAdd = mCkCrd(end:-1:1,:);
                                    mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                                end
                                
                                mTpAct(jCrack,2) = false;
                                mCk2Up(jCrack,2) = true;
                                mTpRdi(jCrack,2) = 0;
                                
                            else
                                
                                if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                    
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCkCrd-i_sgm>3
                                            mCkAdd = mCkCrd(i_sgm+1:end-2,:);
                                        else
                                            mCkAdd = mCkCrd(i_sgm+1:end-1,:);
                                        end
                                    else
                                        % mCkAdd = mCkCrd(i_sgm+1:end,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm+1:end,:);
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,1) % just incremented
                                            if mCkJun(jCrack,1) && i_sgm>2
                                                mCkAdd = mCkCrd(i_sgm+1:-1:3,:);
                                            else
                                                mCkAdd = mCkCrd(i_sgm+1:-1:2,:);
                                            end
                                        else
                                            % mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                            mCkAdd = cCkCrd{jCrack}(j_sgm+1:-1:1,:);
                                        end
                                    end
                                    
                                else
                                    
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>2
                                            mCkAdd = mCkCrd(i_sgm+1:-1:3,:);
                                        else
                                            mCkAdd = mCkCrd(i_sgm+1:-1:2,:);
                                        end
                                    else
                                        % mCkAdd = mCkCrd(i_sgm+1:-1:1,:);
                                        mCkAdd = cCkCrd{jCrack}(j_sgm+1:-1:1,:);
                                    end
                                    if size(mCkAdd,1)<2
                                        % inadequate slave crack deflection direction
                                        % try deflecting jCrack in opposite direction
                                        if mCk2Up(jCrack,2) % just incremented
                                            if mCkJun(jCrack,2) && nCkCrd-i_sgm>3
                                                mCkAdd = mCkCrd(i_sgm+1:end-2,:);
                                            else
                                                mCkAdd = mCkCrd(i_sgm+1:end-1,:);
                                            end
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
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1);
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                            
                            mTpAct(iCrack,2) = false;
                            mTpRdi(iCrack,2) = 0;
                            
                            break
                            
                        end
                        
                    else % xs is between A and B
                        
                        if nCkCrd > 2
                            if d_tip(1)*d_sgm(1)+d_tip(2)*d_sgm(2) > 0
                                
                                if mCk2Up(jCrack,2) % just incremented
                                    if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                        mCkAdd = [x_crs;mCkCrd(i_sgm+1:nCkCrd-2,:)];
                                    else
                                        mCkAdd = [x_crs;mCkCrd(i_sgm+1:end-1,:)];
                                    end
                                else
                                    % mCkAdd = [x_crs;mCkCrd(i_sgm+1:end,:)];
                                    mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm+1:end,:)];
                                end
                                if size(mCkAdd,1)<2
                                    % inadequate slave crack deflection direction
                                    % try deflecting jCrack in opposite direction
                                    if mCk2Up(jCrack,1) % just incremented
                                        if mCkJun(jCrack,1) && i_sgm>2
                                            mCkAdd = [x_crs;mCkCrd(i_sgm:-1:3,:)];
                                        else
                                            mCkAdd = [x_crs;mCkCrd(i_sgm:-1:2,:)];
                                        end
                                    else
                                        % mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                                        mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm:-1:1,:)];
                                    end
                                end
                                
                            else
                                
                                if mCk2Up(jCrack,1) % just incremented
                                    if mCkJun(jCrack,1) && i_sgm>2
                                        mCkAdd = [x_crs;mCkCrd(i_sgm:-1:3,:)];
                                    else
                                        mCkAdd = [x_crs;mCkCrd(i_sgm:-1:2,:)];
                                    end
                                else
                                    % mCkAdd = [x_crs;mCkCrd(i_sgm:-1:1,:)];
                                    mCkAdd = [x_crs;cCkCrd{jCrack}(j_sgm:-1:1,:)];
                                end
                                if size(mCkAdd,1)<2
                                    % inadequate slave crack deflection direction
                                    % try deflecting jCrack in opposite direction
                                    if mCk2Up(jCrack,2) % just incremented
                                        if mCkJun(jCrack,2) && nCkCrd-i_sgm>2
                                            mCkAdd = [x_crs;mCkCrd(i_sgm+1:nCkCrd-2,:)];
                                        else
                                            mCkAdd = [x_crs;mCkCrd(i_sgm+1:end-1,:)];
                                        end
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
                        
                        mNoInc(iCrack,2) = size(mCkAdd,1);
                        mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                        cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                        
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
