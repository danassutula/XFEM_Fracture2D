
%--------------------------------------------------------------------------
% Generate mesh (shape: rectangle)
%--------------------------------------------------------------------------

% dimensions
L = 1.0; % job_domain_L;
H = 1.0; % job_domain_H;

if with_Layered
    
    % e.g. yi = [-H/2,-H/4,H/4,H/2]; (size = n_layers+1)
    yi = [-H/2,-H/4,H/4,H/2]; % job_domain_yi; 
    
    % e.g. yi = [0,2,0]; (size = n_layers)
    Ni = [0,2,0]; % job_domain_Ni;
    
end

% he = 1e-2; % job_mesh_hElmStd;

nH = 50; % job_mesh_nElmStd;
ne = ceil([L/H,1]*nH);
he = H/ne(2);

xlim = [   0,  L];
ylim = [-H/2,H/2];

switch mesh_ElemType
    case 'T3'
        if with_Layered
            
            [mNdCrd,mLNodS,cBCNod,vElPhz] = MeshLayers_t3(xlim,yi,Ni,he);
            if any(any(mLNodS)<1); error('Mesh too coarse for some layers'); end
            
        else
            
            [mNdCrd,mLNodS,cBCNod] = MeshRect_t3(xlim,ylim,he);
            vElPhz = ones(size(mLNodS,1),1);
            
        end
    case 'Q4'
        if with_Layered
            
            [mNdCrd,mLNodS,cBCNod,vElPhz] = MeshLayers_q4(xlim,yi,Ni,he);
            if any(any(mLNodS)<1); error('Mesh too coarse for some layers'); end
            
        else
            
            [mNdCrd,mLNodS,cBCNod] = MeshRect_q4(xlim,ylim,ne(1),ne(2));
            vElPhz = ones(size(mLNodS,1),1);
            
        end
end

% boundary vertix nodes
cBCNod{5} = cBCNod{1}(1);
cBCNod{6} = cBCNod{1}(end);
cBCNod{7} = cBCNod{3}(end);
cBCNod{8} = cBCNod{3}(1);

% BC coordinates (alt.)
cBCCrd = [];

clear xlim ylim yi Ni

%--------------------------------------------------------------------------
