function [mNoInc,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun] = GeoCrk_Intersect_tipdist ...
    (qCk2Up,mCk2Up,cCkCrd,mTpAct,mTpRdi,mCkJun,nSgBln,f_tip,f_xrs)

% tol. for min. dist. to ck. vtx. 
f_vtx = 0.1; % as frac. of 'he'

% max n. bln. sgm. to add
nSgAdd = nSgBln; % 3; inf;

msg0 = 'a crack should not intersect another crack for the second time; ';
msg1 = 'removing tip enrichment; i_crk = %d, i_tip = 1'; msg1=[msg0,msg1];
msg2 = 'removing tip enrichment; i_crk = %d, i_tip = 2'; msg2=[msg0,msg2];

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

for iCk2Ch = qCk2Up'
    
    iCrack = iCk2Ch(1);
    iCkTip = iCk2Ch(2);
    
    if iCkTip == 1 && mTpRdi(iCrack,1)
        
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
            
            % get distance from vTpCrd_ref to other crack mCkCrd
            [r,d_inc,i_sgm] = LevelSet2Poly_abs(vTpCrd_ref,mCkCrd);
            
            % tip segment direction
            d_tip = vTpCrd_ref-mCkCrd_ref(2,:);
            
            if r<tol_xrs && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1))
                    
                    % A ----x--> B
                    
                    x_crs = vTpCrd_ref + d_inc;
                    
                    x_ptA = mCkCrd(i_sgm,:);
                    x_ptB = mCkCrd(i_sgm+1,:);
                    
                    d_sgm = x_ptB - x_ptA;
                    d_xsg = x_crs - x_ptA;
                    
                    s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                    s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                    
                    if i_sgm == 1 && s_xsg < tol_vtx % xs is on tip 1
                        if mCk2Up(iCrack,1) && mTpAct(iCrack,1) % incremented tip
                            
                            if mCk2Up(jCrack,2) % incremented tip
                                if mCkJun(jCrack,2) % just intersected
                                    mCkAdd = mCkCrd(end-2:-1:1,:);
                                else
                                    mCkAdd = mCkCrd(end-1:-1:1,:);
                                end
                            else
                                % mCkAdd = mCkCrd(end:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
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
                            
                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;
                            
                        end
                        
                        mTpAct(iCrack,1) = false;
                        mCk2Up(iCrack,1) = true;
                        mTpRdi(iCrack,1) = 0;
                        
                        break
                        
                    elseif i_sgm+1 == nCkCrd && s_sgm-s_xsg < tol_vtx % xs on tip 2
                        if mCk2Up(iCrack,1) && mTpAct(iCrack,1) % incremented tip
                            
                            if mCk2Up(jCrack,1) % incremented tip
                                if mCkJun(jCrack,1) % just intersected
                                    mCkAdd = mCkCrd(3:end,:);
                                else
                                    mCkAdd = mCkCrd(2:end,:);
                                end
                            else
                                % mCkAdd = mCkCrd;
                                mCkAdd = cCkCrd{jCrack};
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
                            
                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;
                            
                        end
                        
                        mTpAct(iCrack,1) = false;
                        mCk2Up(iCrack,1) = true;
                        mTpRdi(iCrack,1) = 0;
                        
                        break
                            
                    end
                    
                else
                    
                    warning(msg1,iCrack);
                    mTpAct(iCrack,1) = false;
                    mTpRdi(iCrack,1) = 0;
                    
                    break
                    
                end
            end
        end
        
    elseif iCkTip == 2 && mTpRdi(iCrack,2) > 0
        
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
            
            % get distance from vTpCrd_ref to other crack mCkCrd
            [r,d_inc,i_sgm] = LevelSet2Poly_abs(vTpCrd_ref,mCkCrd);
            
            % tip segment direction
            d_tip = vTpCrd_ref-mCkCrd_ref(end-1,:);
            
            if r<tol_xrs && d_inc(1)*d_tip(1)+d_inc(2)*d_tip(2) >= 0
                if isempty(GeoXPoly(mCkCrd,mCkCrd_ref,1)) 
                    
                    % A ----x--> B
                    
                    x_crs = vTpCrd_ref + d_inc;
                    
                    x_ptA = mCkCrd(i_sgm,:);
                    x_ptB = mCkCrd(i_sgm+1,:);
                    
                    d_sgm = x_ptB - x_ptA;
                    d_xsg = x_crs - x_ptA;
                    
                    s_sgm = sqrt(d_sgm(1)^2+d_sgm(2)^2);
                    s_xsg = sqrt(d_xsg(1)^2+d_xsg(2)^2);
                    
                    if i_sgm == 1 && s_xsg < tol_vtx % xs is on tip 1
                        if mCk2Up(iCrack,2) && mTpAct(iCrack,2) % incremented tip
                            
                            if mCk2Up(jCrack,2) % just incremented
                                if mCkJun(jCrack,2) % just intersected
                                    mCkAdd = mCkCrd(1:end-2,:);
                                else
                                    mCkAdd = mCkCrd(1:end-1,:);
                                end
                            else
                                % mCkAdd = mCkCrd;
                                mCkAdd = cCkCrd{jCrack};
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(1:nSgAdd+1,:);
                            end
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1);
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                            
                            mTpAct(jCrack,1) = false;
                            mCk2Up(jCrack,1) = true;
                            mTpRdi(jCrack,1) = 0;
                            
                        end
                        
                        mTpAct(iCrack,2) = false;
                        mCk2Up(iCrack,2) = true;
                        mTpRdi(iCrack,2) = 0;
                        
                        break
                        
                    elseif i_sgm+1 == nCkCrd && s_sgm-s_xsg < tol_vtx % xs on tip 2
                        if mCk2Up(iCrack,2) && mTpAct(iCrack,2) % incremented tip
                            
                            if mCk2Up(jCrack,1) % just incremented
                                if mCkJun(jCrack,1) % just intersected
                                    mCkAdd = mCkCrd(end:-1:3,:);
                                else
                                    mCkAdd = mCkCrd(end:-1:2,:);
                                end
                            else
                                % mCkAdd = mCkCrd(end:-1:1,:);
                                mCkAdd = cCkCrd{jCrack}(end:-1:1,:);
                            end
                            
                            if size(mCkAdd,1) > nSgAdd+1
                                mCkAdd = mCkAdd(1:nSgAdd+1,:);
                            end
                            
                            mNoInc(iCrack,2) = size(mCkAdd,1);
                            mCkJun(iCrack,2) = length(cCkCrd{iCrack})+1;
                            cCkCrd{iCrack} = [cCkCrd{iCrack};mCkAdd];
                            
                            mTpAct(jCrack,2) = false;
                            mCk2Up(jCrack,2) = true;
                            mTpRdi(jCrack,2) = 0;
                            
                        end
                        
                        mTpAct(iCrack,2) = false;
                        mCk2Up(iCrack,2) = true;
                        mTpRdi(iCrack,2) = 0;
                        
                        break
                        
                    end
                else
                    
                    warning(msg2,iCrack);
                    mTpAct(iCrack,2) = false;
                    mTpRdi(iCrack,2) = 0;
                    
                    break
                    
                end
            end
        end
    end
end
