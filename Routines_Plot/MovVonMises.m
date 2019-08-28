
%==========================================================================
% Plot Von Mises (from file)
%==========================================================================

if iStep <= 1
   fprintf('\niStep = %u\n\n',iStep)
   warning('a movie requires that iStep > 1; no movie was produced.')
   return
end

%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------

% movie duration
T = 10; % (sec.)
% movie frame rate
fps = ceil(iStep/T);
% fps = 3;

deformed_scale = 10

% dot size
dotsz_std = 6; % 12
dotsz_enr = 1; % 3
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


%--------------------------------------------------------------------------
% Figure format
%--------------------------------------------------------------------------

if exist('fig_ssM','var') && exist('szfig_ssM','var')
    h = figure(fig_ssM); set(h,'Color','w');
    set(fig_ssM,'OuterPosition',szfig_ssM);
else
    h = figure; set(h,'Color','w');
    scrsz = get(0,'ScreenSize'); scrsz=scrsz([3,4])-scrsz([1,2])+1;
    set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-360,1280,720]); % wide
    % set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-640,1280,1280]); % square
end

%--------------------------------------------------------------------------
% Normalise stresses
%--------------------------------------------------------------------------

% init. stress bounds:
ss_maxmin = 0; % largest minimum
ss_minmax = inf; % smallest maximum

% init. stress bounds
ss_maxmax = 0; % largest maximum
ss_minmin = inf; % smallest minimum

for i_frm = 1:iStep+1
    
    load([path_savedVar,'var_Stress',num2str(i_frm,'_%05d')]);
    vMises = [sStres.vMises_std;sStres.vMises_enr];
    
    ss_i = max(vMises);
    
    if ss_minmax > ss_i
        ss_minmax = ss_i;
    end
    if ss_maxmax < ss_i
        ss_maxmax = ss_i;
    end
    
    ss_i = min(vMises);
    
    if ss_maxmin < ss_i
        ss_maxmin = ss_i;
    end
    if ss_minmin > ss_i
        ss_minmin = ss_i;
    end
    
end

%--------------------------------------------------------------------------
% Make movie
%--------------------------------------------------------------------------

if exist('save_DisplcAll','var') && save_DisplcAll
    
    load([path_savedVar,'var_Displc',num2str(iStep,'_%05d')]);
    mGsDfm = sDisps.mGsCrd_std + sDisps.mGsDsp_std*deformed_scale;
    
    flag_plotDisplacements = 1;
    flag_plotCracks = 0;
    
else
    
    mGsDfm = sStres.mGsCrd_std;
    
    flag_plotDisplacements = 0;
    
    if exist('save_CracksAll','var') && save_CracksAll
        flag_plotCracks = 1;
    else
        flag_plotCracks = 0;
    end
    
end


xmin = min(mGsDfm(:,1));
xmax = max(mGsDfm(:,1));
ymin = min(mGsDfm(:,2));
ymax = max(mGsDfm(:,2));

m = max(xmax-xmin,ymax-ymin) * 0.05;
abox = [xmin-m,xmax+m,ymin-m,ymax+m];


% easily to initialize
% sMovie = getframe(gcf);
fprintf('\nMovie frames:\n')

movieWriteTime = clock;
movieWriteTime = ['_',num2str(movieWriteTime(4)),...
    '_',num2str(movieWriteTime(5))];

sMovie = VideoWriter([path_savedMov,...
    'mov_VonMises',movieWriteTime,'.avi']);

sMovie.FrameRate = fps;
sMovie.Quality = 100;

open(sMovie);

for i_frm = 1:iStep+1
    
    fprintf('step = %i/%i\n',i_frm,iStep);
    load([path_savedVar,'var_Stress',num2str(i_frm,'_%05d')]);
    
    if flag_plotDisplacements
        load([path_savedVar,'var_Displc',num2str(i_frm,'_%05d')]); % sDisps
    else
        sDisps = struct('mGsDsp_std',0,'mGsDsp_enr',0);
    end
    
    clf(h); axis equal; axis(abox); hold on;
    set(gca,'layer','top','box','on');
    
    title({'Fracture process',num2str([i_frm,iStep],'step %i/%i')});
    
    if exist('lengthUnits','var')
        xlabel(['x (',lengthUnits,')'])
        ylabel(['y (',lengthUnits,')'])
    end
    
    set(gca,'FontSize',szfnt);
    
    
    vMises = sStres.vMises_std;
    mGsCrd = sStres.mGsCrd_std;
    mGsDsp = sDisps.mGsDsp_std;
    
    vMises(vMises<ss_maxmin) = ss_maxmin;
    vMises = log10(vMises./ss_maxmax); % scaling
    
    mGsDfm = mGsCrd + mGsDsp*deformed_scale;
    scatter(mGsDfm(:,1),mGsDfm(:,2),dotsz_std,vMises,'o','fill')
    
    
    vMises = sStres.vMises_enr;
    mGsCrd = sStres.mGsCrd_enr;
    mGsDsp = sDisps.mGsDsp_enr;
    
    vMises(vMises<ss_maxmin) = ss_maxmin;
    vMises = log10(vMises./ss_maxmax); % scaling
    
    mGsDfm = mGsCrd + mGsDsp*deformed_scale;
    scatter(mGsDfm(:,1),mGsDfm(:,2),dotsz_enr,vMises,'o','fill')
    
    
    if flag_plotCracks
        
        load([path_savedVar,'var_Crack',num2str(i_frm,'_%05d')]);
        cCkCrd = sCrack.cCkCrd;
        
        if flag_plotCracks && exist('cCkCrd','var')
            for j = 1:length(cCkCrd)
                plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'w','linewidth',1)
            end
        end
        
    end
    
    
    colormap(jet); % bone, hot, cool
    c = colorbar('location','EastOutside');
    c.Label.String = 'log_{10} (\sigma_{vms}/\sigma_{vmsMax})';
    
    % sMovie(ii) = getframe(gcf);
    currentMovieFrame = getframe(gcf);
    writeVideo(sMovie,currentMovieFrame.cdata);
    
end

% movie2avi(sMovie(1:ii),[path_savedMov,'mov_VonMises',tmp],...
%    'fps',fps) % ,'compression','none'

close(sMovie)

clear vMises mGsCrd mGsDsp mGsDfm
clear ss_maxmax ss_maxmin ss_minmax ss_minmin
clear ss_i q_std q_enr abox m

clear flag_plotDisplacements
clear flag_plotCracks

fprintf('\nMovie done.\n')

close(h)
