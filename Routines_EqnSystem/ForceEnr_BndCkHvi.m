function Fg = ForceEnr_BndCkHvi(nGlDof,nCrack,cCkCrd,N_gam,W,P)

% P(i_crk,:) = [ p , \tau ]

%--------------------------------------------------------------------------
% Call Global Variables
%--------------------------------------------------------------------------

global mNdCrd
global mLNodS

global cLEnDt
global cLNodE
global mNDofE

global cHvStd_elm
global cHvStd_sgm
global cHvBln_elm
global cHvBln_sgm

global with_MapTyp
global tol_abs

%--------------------------------------------------------------------------

nLNodS = length(mLNodS(1,:));

n_elD = nLNodS*2;
n_spG = nGlDof*2; % approx.
n_gsp = length(W);

ig = zeros(n_spG,1);
Fg = zeros(n_spG,1);

p_ndS = 1:nLNodS;
p_map = ones(nLNodS,1);

q = 1:n_elD;

for i_crk = 1:nCrack
    p_crk = P(i_crk,:);
    
    if any(p_crk)
        
        x_crk = cCkCrd{i_crk};
        
        % std. & bln. (no bln. ramp)
        vElEnr = [cHvBln_elm{i_crk,1};cHvStd_elm{i_crk};cHvBln_elm{i_crk,2}];
        mCkSgm = [cHvBln_sgm{i_crk,1};cHvStd_sgm{i_crk};cHvBln_sgm{i_crk,2}];
        
        for i_elE = 1:length(vElEnr)
            
            i_elm = vElEnr(i_elE);
            x_elm = mNdCrd(mLNodS(i_elm,:),:);
            
            % element enr. data
            mLEnDt = cLEnDt{i_elm};
            
            i_enr = find(mLEnDt(1,:) == i_crk);
            i_enr = i_enr(mLEnDt(2,i_enr) == 0);
            
            % get element enriched nodes and DOF 
            p_ndE = p_ndS + sum(mLEnDt(3,1:i_enr-1));
            p_dof = mNDofE(:,cLNodE{i_elm}(p_ndE));
            
            for i_sgm = mCkSgm(i_elE,1):mCkSgm(i_elE,2)
                
                x_sgm = x_crk([i_sgm,i_sgm+1],:);
                [n_crs,x_crs] = GeoXElem(x_sgm,x_elm);
                
                % to integrate, must have a length and must have a jump 
                if n_crs == 2 && ~all(LevelSet2Line(x_elm,x_sgm)>-tol_abs)
                    
                    n = x_crs(2,:) - x_crs(1,:);
                    s = sqrt(n(1)^2+n(2)^2);
                    n = s\n; % detJ = 0.5*s;
                    
                    f = [n(1)*p_crk(2) - n(2)*p_crk(1);
                         n(1)*p_crk(1) + n(2)*p_crk(2)];
                    
                    x_gsp = N_gam*x_crs;
                    
                    if with_MapTyp == 1
                        % local mapping
                        X_gsp = N_gam*Gauss_Glb2Lcl(x_crs,x_elm);
                    else
                        % global mapping
                        X_gsp = Gauss_Glb2Lcl(x_gsp,x_elm);
                    end
                    
                    H = EnrJmp_Hvi(x_gsp); % n_gsp x 1
                    N = ShapesStd_omg(n_gsp,X_gsp);
                    F = f*sum(N.*W(:,p_map)*(H(1)*0.5*s),1);
                    
                    ig(q) = p_dof(:);
                    Fg(q) = F(:);
                    
                    q = q + n_elD;
                    
                end
            end
        end
    end
end

q=1:q(1)-1;
ig = ig(q);
Fg = Fg(q);

Fg = sparse(ig,ones(length(ig),1),Fg,nGlDof,1);

end
