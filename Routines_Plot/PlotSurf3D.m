function PlotSurf3D(coord,lnods)

h = gcf; clf(h); set(h,'Color','w')

[n_elm,n_nod] = size(lnods);

x = zeros(n_nod,n_elm); 
y = zeros(n_nod,n_elm);
z = zeros(n_nod,n_elm);

for i = 1:n_elm
    x(:,i) = coord(lnods(i,:),1);
    y(:,i) = coord(lnods(i,:),2);
    z(:,i) = coord(lnods(i,:),3);
end

patch(x,y,z,z,'EdgeColor','none');

axis equal
box on

xmin = min(coord(:,1));
xmax = max(coord(:,1));
ymin = min(coord(:,2));
ymax = max(coord(:,2));
zmin = min(coord(:,3));
zmax = max(coord(:,3));

m = max([xmax-xmin,ymax-ymin,zmax-zmin])*0.05; % 5% margins
axis([xmin-m,xmax+m,ymin-m,ymax+m,zmin-m,zmax+m])
set(gca,'layer','top');

%==========================================================================

end

