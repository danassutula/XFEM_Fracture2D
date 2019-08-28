
%==========================================================================
% Plot Cracks (from file)
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
T = 6; % (sec.)
% movie frame rate
fps = round(iStep/T);
% text font size
szfnt = 22;

% limit frame rate
if fps>12; fps=12;
    T = iStep/fps;
end

%--------------------------------------------------------------------------
% Figure format
%--------------------------------------------------------------------------

if exist('fig_ckM','var') && exist('szfig_ckM','var')
    h = figure(fig_ckM); set(fig_ckM,'OuterPosition',szfig_ckM);
else
    h = figure; set(h,'Color','w');
    scrsz = get(0,'ScreenSize'); scrsz=scrsz([3,4])-scrsz([1,2])+1;
    set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-360,1280,720]); % wide
    % set(h,'OuterPosition',[scrsz(1)/2-640,scrsz(2)/2-640,1280,1280]); % square
end

%--------------------------------------------------------------------------
% Make movie
%--------------------------------------------------------------------------

% PlotDomain;

fprintf('\nMovie frames:\n')

movieWriteTime = clock;
movieWriteTime = ['_',num2str(movieWriteTime(4)),...
    '_',num2str(movieWriteTime(5))];

sMovie = VideoWriter([path_savedMov,...
    'mov_CrackGrowth',movieWriteTime,'.avi']);

sMovie.FrameRate = fps;
sMovie.Quality = 100;

open(sMovie);

for i_frm = 1:iStep+1
    
    fprintf('step = %i/%i\n',i_frm,iStep);
    load([path_savedVar,'var_Crack',num2str(i_frm,'_%05d')]);
    
    cCkCrd = sCrack.cCkCrd;
    mTpAct = sCrack.mTpAct;
    
    PlotDomain
    PlotCracks
    
    title({'Fracture process',num2str([i_frm,iStep],'step %i/%i')});
    
    set(gca,'FontSize',szfnt);
    
    % sMovie(i_frm) = getframe(gcf);
    currentMovieFrame = getframe(gcf);
    writeVideo(sMovie,currentMovieFrame.cdata);
    
end

close(sMovie)

clear sMovie

fprintf('\nMovie done.\n')

close(h)

%--------------------------------------------------------------------------
