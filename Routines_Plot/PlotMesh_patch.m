function PlotMesh_patch(elems)

global mNdCrd
global mLNodS

for i = elems(:)'
    patch(mNdCrd(mLNodS(i,:),1),mNdCrd(mLNodS(i,:),2),...
        [0,0.75,0.75],'facealpha',0.5,'edgecolor','b');
end