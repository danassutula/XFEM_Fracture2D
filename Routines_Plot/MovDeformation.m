
%==========================================================================
% MovDeformation
%==========================================================================

if iStep <= 1
   fprintf('\niStep = %u\n\n',iStep)
   warning('a movie requires that iStep > 1; no movie was produced.')
   return
end

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

FigFormat 

% movie duration
T = 20; % (sec.)
% movie frame rate
fps = ceil(iStep/T);
% fps = 3;

% dot size
dotsz_std = 3;
dotsz_enr = 1;

% text font size
szfnt = 22;

if ~exist('deformed_scale','var')
    deformed_scale = 10;
end

if ~exist('path_savedVar','var')
    path_savedVar = [cd,'\'];
end

% limit frame rate
if fps>12; fps=12;
    T = iStep/fps;
end

deformed_scale = 5

%--------------------------------------------------------------------------
% Figure format
%--------------------------------------------------------------------------

if exist('fig_dfM','var') && exist('szfig_dfM','var')
    h = figure(fig_dfM); set(h,'Color','w');
    set(fig_dfM,'OuterPosition',szfig_dfM);
else
    h = figure; set(h,'Color','w');
    scrsz = get(0,'ScreenSize'); scrsz=scrsz([3,4])-scrsz([1,2])+1;
    set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-360,1280,720]); % wide
    % set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-640,1280,1280]); % square
end

%--------------------------------------------------------------------------
% Make movie
%--------------------------------------------------------------------------

load([path_savedVar,'var_Displc',num2str(1,'_%05d')]);

mGsCrd = sDisps.mGsCrd_std;
mGsDsp = sDisps.mGsDsp_std;
mGsDfm = mGsCrd + mGsDsp*deformed_scale;

% % ALT-1
% xmin = min(mGsCrd(:,1));
% xmax = max(mGsCrd(:,1));
% ymin = min(mGsCrd(:,2));
% ymax = max(mGsCrd(:,2));

% ALT-2
xmin = min(mGsDfm(:,1));
xmax = max(mGsDfm(:,1));
ymin = min(mGsDfm(:,2));
ymax = max(mGsDfm(:,2));

m = max(xmax-xmin,ymax-ymin) * 0.05;
abox = [xmin-m,xmax+m,ymin-m,ymax+m];

% easily to initialize
sMovie = getframe(gcf);
fprintf('\nMovie frames:\n')

for i_frm = 1:iStep+1
    
    fprintf('step = %i/%i\n',i_frm,iStep);
    load([path_savedVar,'var_Displc',num2str(i_frm,'_%05d')]);
    
    
    clf(h); axis equal; axis(abox); hold on;
    set(gca,'layer','top','box','on');
    
    title({'Fracture process',num2str([i_frm,iStep],'step %i/%i')});
    
    if exist('lengthUnits','var')
        xlabel(['x (',lengthUnits,')'])
        ylabel(['y (',lengthUnits,')'])
    end
    
    set(gca,'FontSize',szfnt);
    
    
    mGsCrd = sDisps.mGsCrd_std;
    mGsDsp = sDisps.mGsDsp_std;
    mGsDfm = mGsCrd + mGsDsp*deformed_scale;
    scatter(mGsDfm(:,1),mGsDfm(:,2),dotsz_std,'k','o','fill')
    
    pause(1/3)
    
    mGsCrd = sDisps.mGsCrd_enr;
    mGsDsp = sDisps.mGsDsp_enr;
    mGsDfm = mGsCrd + mGsDsp*deformed_scale;
    scatter(mGsDfm(:,1),mGsDfm(:,2),dotsz_enr,'r','o','fill')
    
    pause(1/3)
    
    
    sMovie(i_frm) = getframe(gcf);
    
end

tmp = clock; tmp = ['_',num2str(tmp(4)),'_',num2str(tmp(5))];

movie2avi(sMovie(1:i_frm),[path_savedMov,'mov_Deformation',tmp],...
    'fps',fps) % ,'compression','none'

% clear sMovie 
clear mGsCrd mGsDsp q_std q_enr abox m

close(h)
