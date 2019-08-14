
%==========================================================================
% Plot Enriched Elements
%==========================================================================

gcf; clf;
PlotMesh;
hold on; 

% init. patch coord.
x = zeros(nLNodS,nElEnr);
y = zeros(nLNodS,nElEnr);

c_hvi = [];
c_brn = [];

% patch colors
c_hvi(:,1) = linspace(0,7/11,7) * 1;
c_hvi(:,2) = linspace(0,7/11,7) * 1;
c_hvi(:,3) = linspace(0,7/11,7) * 1;
c_hvi = c_hvi(3:end,:);

c_brn(:,1) = linspace(7/11,1.0,5) * 1;
c_brn(:,2) = linspace(7/11,1.0,5) * 1;
c_brn(:,3) = linspace(7/11,1.0,5) * 1;
c_brn = c_brn(2:end-1,:);

%--------------------------------------------------------------------------
% Branch
%--------------------------------------------------------------------------

n = 0;

for i = 1:nCrack
    for j = 1:2
        for k = cBrStd_eTp{i,j}'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(k,:),1);
            y(:,n) = mNdCrd(mLNodS(k,:),2);
            
        end
    end
end

h1 = patch(x(:,1:n),y(:,1:n),c_brn(1,:),'edgecolor','k','linestyle','-','linewidth',1.75);

n = 0;

for i = 1:nCrack
    for j = 1:2
        for k = cBrStd_eFl{i,j}'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(k,:),1);
            y(:,n) = mNdCrd(mLNodS(k,:),2);
            
        end
    end
end

h2 = patch(x(:,1:n),y(:,1:n),c_brn(2,:));

n = 0;

for i = 1:nCrack
    for j = 1:2
        for k = cBrBln_eFl{i,j}'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(k,:),1);
            y(:,n) = mNdCrd(mLNodS(k,:),2);
            
        end
    end
end

h3 = patch(x(:,1:n),y(:,1:n),c_brn(3,:),'edgecolor','k','linestyle','--','linewidth',1.5);

%--------------------------------------------------------------------------
% Heaviside
%--------------------------------------------------------------------------

n = 0;

for i = 1:nCrack
    for j = [cHvStd_elm{i}]'
        
        n = n + 1;
        
        x(:,n) = mNdCrd(mLNodS(j,:),1);
        y(:,n) = mNdCrd(mLNodS(j,:),2);
        
    end
end

h4 = patch(x(:,1:n),y(:,1:n),c_hvi(1,:));

n = 0;

for i = 1:nCrack
    for k = 1:2
        for j = cHvStd_elm{i}(ismember(cHvStd_elm{i},cBrBln_eSp{i,k}))'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(j,:),1);
            y(:,n) = mNdCrd(mLNodS(j,:),2);
            
        end
    end
end

h5 = patch(x(:,1:n),y(:,1:n),c_hvi(2,:));

n = 0;

for i = 1:nCrack
    for k = 1:2
        for j = cHvBln_elm{i,k}(ismember(cHvBln_elm{i,k},cBrBln_eSp{i,k}))'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(j,:),1);
            y(:,n) = mNdCrd(mLNodS(j,:),2);
            
        end
    end
end

h6 = patch(x(:,1:n),y(:,1:n),c_hvi(3,:));

n = 0;

for i = 1:nCrack
    for k = 1:2
        for j = cHvStd_elm{i}(ismember(cHvStd_elm{i},cBrStd_eSp{i,k}))'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(j,:),1);
            y(:,n) = mNdCrd(mLNodS(j,:),2);
            
        end
    end
end

h7 = patch(x(:,1:n),y(:,1:n),c_hvi(4,:));

n = 0;

for i = 1:nCrack
    for k = 1:2
        for j = cHvBln_elm{i,k}(ismember(cHvBln_elm{i,k},cBrStd_eSp{i,k}))'
            
            n = n + 1;
            
            x(:,n) = mNdCrd(mLNodS(j,:),1);
            y(:,n) = mNdCrd(mLNodS(j,:),2);
            
        end
    end
end

h8 = patch(x(:,1:n),y(:,1:n),c_hvi(5,:));

%--------------------------------------------------------------------------
% Label all entities (if exist)
%--------------------------------------------------------------------------

h={h1,h2,h3,h4,h5,h6,h7,h8};
q_keep=true(1,8);
for i = 1:8
    if isempty(h{i})
        q_keep(i)= ~1;
    end
end

h_labels={...
    'tip element (branch enrichment)',...
    'non-split el. (branch standard enr.)',...
    'non-split el. (branch blending enr.)',...
    'split el. (Heaviside std. enrichment)',...
    'split el. (Heaviside std. & branch bln.)',...
    'split el. (Heaviside bln. & branch bln.)',...
    'split el. (Heaviside std. & branch std.)',...
    'split el. (Heaviside bln. & branch std.)'};
h_labels = h_labels(q_keep);

hl = legend([h1,h2,h3,h4,h5,h6,h7,h8],h_labels,'location','SouthEast');

%--------------------------------------------------------------------------
% Plot Cracks
%--------------------------------------------------------------------------

for iCrack = 1:nCrack
    scatter(cCkCrd{iCrack}(:,1),cCkCrd{iCrack}(:,2),10,'ow','fill')
    plot(cCkCrd{iCrack}(:,1),cCkCrd{iCrack}(:,2),'w','linewidth',1)
end

%--------------------------------------------------------------------------
% Format figure (for LaTeX)
%--------------------------------------------------------------------------

% title('Enriched element classification')
box off
% axis off

FigFormat
FigResize(szfig_big,szfnt_big)
set(hl,'fontsize',17)

SavePDF('enriched_element_classification',gcf)

% clear var.'s
clear x y n c_brn c_hvi hl

%--------------------------------------------------------------------------
