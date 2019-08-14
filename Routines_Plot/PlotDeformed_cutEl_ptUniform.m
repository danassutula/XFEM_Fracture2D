
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

h = gcf; hold on;
set(h,'Color','w');

if ~exist('deformed_scale','var')
    deformed_scale = 10; end
if ~exist('deformed_dotsz','var')
    deformed_dotsz = 3; end

deformed_scale
deformed_dotsz

color_gpt = [1.0,0,0];

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% TEMP.
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

switch mesh_ElemType
    case 'T3'
        
        X_ = linspace(0,1,10);
        X_ = repmat(X_,length(X_),1);
        X_ = [X_(:),reshape(X_',[],1)];
        X_(sum(X_,2)>1,:) = [];
        
    case 'Q4'
        
        X_ = linspace(-1,1,10);
        X_ = repmat(X_,length(X_),1);
        X_ = [X_(:),reshape(X_',[],1)];
        
    otherwise
        error('mesh_ElemType = ?')
end

N_s = ShapesStd_omg(size(X_,1),X_);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



%--------------------------------------------------------------------------
% Plotting
%--------------------------------------------------------------------------

elems_cut = unique([cell2mat(cBrStd_eTp(:)); ...
    cell2mat(cHvBln_elm(:));cell2mat(cHvStd_elm(:))]);
% elems_nct = setdiff([vElStd(:);vElEnr(:)],elems_cut);

if ~isempty(vElEnr) && 1
    mGsDsp = cell(nElEnr,1);
    
    for ii = 1:length(elems_cut) % vElEnr
        jj = elems_cut(ii); % vElEnr(ii);
        
        vLNodS = mLNodS(jj,:);
        vLNodE = cLNodE{jj};
        
        mElCrd = mNdCrd(vLNodS,:);
        mGsCrd = N_s*mElCrd;
        
        mLDspS = mNDspS(:,vLNodS)*deformed_scale;
        mLDspE = mNDspE(:,vLNodE)*deformed_scale;
        
        mGsDsp{ii} = N_s*(mElCrd+mLDspS');
        
        N_e = [];
            
        for kk = 1:size(cLEnDt{jj},2)
            
            i_crk = cLEnDt{jj}(1,kk);
            switch cLEnDt{jj}(2,kk)
                case 0
                    
                    [H_gsp,H_nod] = EnrFun_Hvi(mElCrd,mGsCrd,cCkCrd{i_crk});
                    N_e = [N_e,N_s.*(repmat(H_gsp,1,size(N_s,2))-repmat(H_nod,size(N_s,1),1))];
                    
                case 1
                    
                    v_tip = diff(cCkCrd{i_crk}([2,1],:));
                    b_tip = atan2(v_tip(2),v_tip(1));
                    
                    [B_gsp,~,B_nod] = EnrFun_Brn(mElCrd,mGsCrd,nEnFun_brn,cCkCrd{i_crk}(1,:),b_tip);
                    
                    if any(jj == cBrBln_eSp{i_crk,1})
                        jj_rmp = jj==cBrBln_eSp{i_crk,1};
                        
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,repmat(N_s*cBrBln_rSp{i_crk,1}(jj_rmp,:)',1,size(N_s,2)).*N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    elseif any(jj == cBrBln_eFl{i_crk,1})
                        jj_rmp = jj==cBrBln_eFl{i_crk,1};
                        
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,repmat(N_s*cBrBln_rFl{i_crk,1}(jj_rmp,:)',1,size(N_s,2)).*N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    else
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    end
                    
                case 2
                    
                    v_tip = diff(cCkCrd{i_crk}([end-1,end],:));
                    b_tip = atan2(v_tip(2),v_tip(1));
                    
                    [B_gsp,~,B_nod] = EnrFun_Brn(mElCrd,mGsCrd,nEnFun_brn,cCkCrd{i_crk}(end,:),b_tip);
                    
                    if any(jj == cBrBln_eSp{i_crk,2})
                        jj_rmp = jj==cBrBln_eSp{i_crk,2};
                        
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,repmat(N_s*cBrBln_rSp{i_crk,2}(jj_rmp,:)',1,size(N_s,2)).*N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    elseif any(jj == cBrBln_eFl{i_crk,2})
                        jj_rmp = jj==cBrBln_eFl{i_crk,2};
                        
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,repmat(N_s*cBrBln_rFl{i_crk,2}(jj_rmp,:)',1,size(N_s,2)).*N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    else
                        for i_enr = 1:nEnFun_brn
                            N_e = [N_e,N_s ...
                                .*(repmat(B_gsp(:,i_enr),1,size(N_s,2))-repmat(B_nod(i_enr,:),size(N_s,1),1))];
                        end
                    end
                otherwise
                    
                    warning('unrecognised enrichment key')
                    cLEnDt{jj}(2,kk)
                    keyabord
                    
            end
        end
        
        mGsDsp{ii} = mGsDsp{ii} + N_e*mLDspE';
        
    end
    
    mGsDsp = cell2mat(mGsDsp);
    
    if exist('plotBox_deformed','var')
        q = mGsDsp(:,1) > plotBox_deformed(1,1) & mGsDsp(:,1) < plotBox_deformed(2,1) & ...
            mGsDsp(:,2) > plotBox_deformed(1,2) & mGsDsp(:,2) < plotBox_deformed(2,2) ;
        mGsDsp = mGsDsp(q,:);
    end
    
    scatter(mGsDsp(:,1),mGsDsp(:,2),deformed_dotsz,color_gpt,'o','full') % mGsCrd_dsp(:,1),
    
end

%--------------------------------------------------------------------------
% Format axis
%--------------------------------------------------------------------------

axis equal

if exist('lengthUnits','var')
    xlabel(['x (',lengthUnits,')'])
    ylabel(['y (',lengthUnits,')'])
else
    xlabel('x')
    ylabel('y')
end

% clear var.
clear mGsCrd mGsDsp mGsDfm

%--------------------------------------------------------------------------
