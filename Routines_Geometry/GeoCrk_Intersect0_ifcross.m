function [nCrack,cCkCrd,mTpAct,mTpRdi,mCkJun,mCkNew] = ...
    GeoCrk_Intersect0_ifcross(cCkCrd,mTpAct,mTpRdi,f_rmv)

%--------------------------------------------------------------------------
% PART I: get crack intersections
%--------------------------------------------------------------------------

fos = 3;                    	% assume x10 intersections per crack

nCrack = length(cCkCrd);
nCkFos = nCrack*fos;

mXsCrd = zeros(nCkFos,2);       % minor - major crack intersections
mMinor = zeros(nCkFos,2);       % minor cracks (are cut) and segments
mMajor = zeros(nCkFos,2);       % major cracks (not cut) and segments

k = 0;                          % count intersections

for i_min = nCrack:-1:2
    p_min = cCkCrd{i_min};
    for j_min = size(p_min,1)-1:-1:1
        d_min = p_min([j_min,j_min+1],:);
        
        for i_mjr = 1:i_min-1
            p_mjr = cCkCrd{i_mjr};
            for j_mjr = 1:size(p_mjr,1)-1
                d_mjr = p_mjr([j_mjr,j_mjr+1],:);
                
                x = GeoXLine(d_min,d_mjr);
                
                if ~isempty(x)
                    
                    k = k + 1;
                    
                    mXsCrd(k,1) = x(1);
                    mXsCrd(k,2) = x(2);
                    
                    mMinor(k,1) = i_min;
                    mMinor(k,2) = j_min;
                    
                    mMajor(k,1) = i_mjr;
                    mMajor(k,2) = j_mjr;
                    
                end
            end
        end
    end
end

j = 1:k; nCkXrs = k;

mXsCrd = mXsCrd(j,:);
mMinor = mMinor(j,:); % (:,1) - minor cracks, (:,2) - crack segments
mMajor = mMajor(j,:); % (:,1) - major cracks, (:,2) - crack segments

if nCkXrs == 0
    mCkJun = zeros(nCrack,2);
    mCkNew = [];
    return
end

%--------------------------------------------------------------------------
% PART II: modify cracks
%--------------------------------------------------------------------------

iCkNew = nCrack + 1;
nCrack = nCrack + nCkXrs;

% junctions and tracking
mCkJun = zeros(nCrack,2);
mCkNew = zeros(nCkXrs,2);

% resize for new cracks
cCkCrd{nCrack,1} = [];
mTpAct(nCrack,2) = 0;
mTpRdi(nCrack,2) = 0;

% find unique minor cracks
[~,jCkXrs] = unique(mMinor(:,1),'first');
jCkXrs = sort(jCkXrs,'ascend');
% n. times minor cracks are cut
vNoCkX = diff([jCkXrs;nCkXrs+1],1);

for iIdCkX = 1:length(jCkXrs);                                              % loop over each minor crack
    
    uNoCkX      = vNoCkX(iIdCkX);                                           % n. times minor crk. is cut
    iCkXrs      = jCkXrs(iIdCkX);                                           % index to unq. minor crack
    jCkXrs_sb1  = iCkXrs:iCkXrs+uNoCkX-1;                                   % index for minor crack segment(s)
    
    uCkMin      = mMinor(iCkXrs,1);                                         % minor crack
    vSgMin      = mMinor(jCkXrs_sb1,2);                                     % minor unique segment(s)
    
    [~,jSgXrs]  = unique(vSgMin,'first');                                   % need to be sequential; use 'first'
    jSgXrs      = sort(jSgXrs);                                             % index to unq. segments
    
    vNoSgX      = diff([jSgXrs;uNoCkX+1],1);                                % n. intersections per segment
    
    for iIdSgX = 1:length(jSgXrs)                                           % loop through segment intersections
        
        uNoSgX      = vNoSgX(iIdSgX);                                       % n. intersections for segment
        iSgXrs      = jSgXrs(iIdSgX);                                       % index to unique segment
        uSgMin      = vSgMin(iSgXrs);                                       % corresponding crack segment
        
        iCkXrs      = jCkXrs_sb1(iSgXrs);                                   % get index to minor crack
        jCkXrs_sb2  = iCkXrs:iCkXrs+uNoSgX-1;                               % get index to all intersections
        
        x_ref   = cCkCrd{uCkMin}(uSgMin,:);
        d_xrs   = zeros(uNoSgX,2); k = 1;
        
        for i = jCkXrs_sb2
            d_xrs(k,:) = x_ref-mXsCrd(i,:); k = k+1;
        end
        
        [~,k]  = sort(d_xrs(:,1).^2+d_xrs(:,2).^2,'descend');
        jCkXrs_sb2  = jCkXrs_sb2(k);
        
        for iCkXrs = jCkXrs_sb2
            
            vXsCrd = mXsCrd(iCkXrs,:);
            uCkMjr = mMajor(iCkXrs,1);
            uSgMjr = mMajor(iCkXrs,2);
            
            x_ref = cCkCrd{uCkMjr}(uSgMjr+1,:);
            d_ref = x_ref-cCkCrd{uCkMjr}(uSgMjr,:);         % o--->o d_ref
            d_xrs = x_ref-vXsCrd;                           %   x->o d_xrs
            
            if d_xrs(1)^2+d_xrs(2)^2 > 0.25*(d_ref(1)^2+d_ref(2)^2)
                
                i_min = uSgMjr+1:size(cCkCrd{uCkMjr},1);    % o--o/o---->o
                i_new = i_min(end:-1:1);                    % o--o/o<----o
                
            else
                
                i_new = 1:uSgMjr;                           % o---->o/o--o
                i_min = i_new(end:-1:1);                    % o<----o/o--o
                
            end
            
            cCkCrd{iCkNew} = [cCkCrd{uCkMjr}(i_new,:);...
                vXsCrd;cCkCrd{uCkMin}(uSgMin+1:end,:)];     % create new
            
            cCkCrd{uCkMin} = [cCkCrd{uCkMin}(1:uSgMin,:);...
                vXsCrd;cCkCrd{uCkMjr}(i_min,:)];            % update old
            
            if mCkJun(uCkMin,2)
                
                mCkJun(iCkNew,2) = size(cCkCrd{iCkNew},1)-mCkJun(uCkMin,2); % update ( mCkJun(uCkMin,2) is relative to end )
                
            else
                
                mTpRdi(iCkNew,2) = mTpRdi(uCkMin,2);
                mTpRdi(uCkMin,2) = 0;
                
                mTpAct(iCkNew,2) = mTpAct(uCkMin,2);
                mTpAct(uCkMin,2) = 0;
                
            end
            
            mCkJun(iCkNew,1) = length(i_new);
            mCkJun(uCkMin,2) = length(i_min);                               % relative to end (tmp., to be corrected)
            
            mCkNew(iCkXrs,1) = uCkMin;
            mCkNew(iCkXrs,2) = iCkNew;
            
            iCkNew = iCkNew + 1;
            
        end
    end
end

for i = mMinor(jCkXrs,1)'
    mCkJun(i,2) = size(cCkCrd{i},1)-mCkJun(i,2);                           % correcting minor crack blending segments
end

%--------------------------------------------------------------------------
% PART III: remove short intersecting cracks (cannot be handled numerically)
%--------------------------------------------------------------------------

% !!! DOES NOT CONSIDER MASTER CRACK

% if ~all(mMinor(:,1)==mCkNew(:,1))
%     error('Does not correspond!!! (Your gut was right?)')
% else
%     warning('YOU GUT IS STILL RIGHT.')
% end
% 
% j = true(nCrack,1);
% 
% for i = find(mTpRdi(:,1) & mCkJun(:,2)==2)'
%     
%     k = mMajor(i==mMinor(:,1)|i==mCkNew(:,2),1);
%     s = LevelSet2Poly_abs(cCkCrd{i}(1,:),cCkCrd{k});
%     
%     if s<mTpRdi(i,1)*f_rmv; j(i)=false;
%         warning('tip too close to another crack; removing crack i_crk = %d',i)
%     end
%     
% end
% 
% for i = find(mTpRdi(:,2) & mCkJun(:,1)~=0)'
%     if size(cCkCrd{i},1)-mCkJun(i,1) == 2
%         
%         k = mMajor(i==mMinor(:,1)|i==mCkNew(:,2),1);
%         s = LevelSet2Poly_abs(cCkCrd{i}(end,:),cCkCrd{k});
%         
%         if s<mTpRdi(i,2)*f_rmv; j(i)=false;
%             warning('tip too close to another crack; removing crack i_crk = %d',i)
%         end
%         
%     end
% end
% 
% cCkCrd = cCkCrd(j);
% nCrack = sum(j);
% 
% mCkJun = mCkJun(j,:);
% mTpAct = mTpAct(j,:);
% mTpRdi = mTpRdi(j,:);

%--------------------------------------------------------------------------
% PART IV: remove short cracks (cannot be handled numerically)
%--------------------------------------------------------------------------

% !!! DOES NOT CONSIDER MASTER CRACK

% j = true(nCrack,1);
% 
% for i = find(mCkJun(:,2)==0)'
%     if size(cCkCrd{i},1)-mCkJun(i,1) == 2
%         % crack tip x_tip = cCkCrd{i}(end,:)
%         x = cCkCrd{i}(end,:)-cCkCrd{i}(end-1,:);
%         if x(1)^2+x(2)^2 < 2*mTpRdi(i,2)^2; j(i)=false;
%             warning('mesh too coarse; i_crk = %d; removing crack',i)
%         end
%     end
% end
% 
% for i = find(mCkJun(:,1)==0)'
%     if mCkJun(i,2) == 2
%         % crack tip x_tip = cCkCrd{i}(1,:)
%         x = cCkCrd{i}(1,:)-cCkCrd{i}(2,:);
%         if x(1)^2+x(2)^2 < 2*mTpRdi(i,1)^2; j(i)=false;
%             warning('mesh too coarse; i_crk = %d; removing crack',i)
%         end
%     end
% end
% 
% cCkCrd = cCkCrd(j);
% nCrack = sum(j);
% 
% mCkJun = mCkJun(j,:);
% mTpAct = mTpAct(j,:);
% mTpRdi = mTpRdi(j,:);

%--------------------------------------------------------------------------
% PART IV: keep track of intersected cracks
%--------------------------------------------------------------------------

% mCkNew
% mMajor
% mMinor

%--------------------------------------------------------------------------

end