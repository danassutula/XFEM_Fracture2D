
h = gcf; clf(h);
set(h,'Color','w');
hold on;

if ~exist('deformed_scale','var')
    deformed_scale = 10;
end
if ~exist('deformed_dotsz','var')
    deformed_dotsz = 3;
end

if 1
    
    PlotDeformed_nctEl
    PlotDeformed_cutEl_ptGauss 
    
else
    
    PlotDeformed_nctEl
    PlotDeformed_cutEl_ptUniform
    
end

mNdDfm = mNdCrd+mNDspS'*deformed_scale;

xmin = min(mNdDfm(:,1));
ymin = min(mNdDfm(:,2));
xmax = max(mNdDfm(:,1));
ymax = max(mNdDfm(:,2));

m = max(xmax-xmin,ymax-ymin)*0.05;
axis([xmin-m,xmax+m,ymin-m,ymax+m])
set(gca,'layer','top','box','on')

xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

title(sprintf('Deformation (x%i)',deformed_scale))

figure(gcf)
