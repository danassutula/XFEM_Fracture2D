
%==========================================================================
% Plot Cracks (from file)
%==========================================================================
%%
close all

% PlotDomain;

% Figure formating
run([cd,'\other_tools\FigFormat.m'])

path0 = '\Results\waferSplitting_mechanical_compareEnergyMin';
open([cd,path0,'\domain.fig']); axis normal

hold on; pause(0.1); fps = 24;
FigResize(szfig_720p,szfnt_big)

title('Fracture process of \itSi\rm-wafer splitting');



path1 = [cd,path0,'\explicit\var\'];
path2 = [cd,path0,'\implicit\var\'];
pathS = [cd,'\']; % saving path

iStep = 557; % (628)
k     = 0;

speed = -0.0015;

n_zm1 = fps * 2;
n_zm2 = fps * 2;
n_hol = fps * 0.5;

xlim_nml = get(gca,'xlim');
ylim_nml = get(gca,'ylim');

xlim_zm1 = [ 0.336143858128238   0.509701663746556];
ylim_zm1 = [ -0.162793095101941   0.160159855971262] * 1e-3;

xlim_zm2 = [0.448919513709791   0.471686687061846];
ylim_zm2 = [-0.352758134888063   0.352241135383917] * 1e-3;



% load ck. coord.
load([path1,'var_Crack',num2str(0,'_%05d')]);
% cCkCrd = sCrack.cCkCrd;

for j = 1:length(cCkCrd)
    a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',1);
end
% legend(a0,'pre-existing flaws')
error
x_min = linspace(xlim_nml(1),xlim_zm1(1),n_zm1);
x_max = linspace(xlim_nml(2),xlim_zm1(2),n_zm1);

y_min = linspace(ylim_nml(1),ylim_zm1(1),n_zm1);
y_max = linspace(ylim_nml(2),ylim_zm1(2),n_zm1);

for i = 1:n_zm1; k=k+1;
    set(gca,'xlim',[x_min(i),x_max(i)])
    set(gca,'ylim',[y_min(i),y_max(i)])
    pause(0.1); sMovie(k) = getframe(gcf);
end

for i = 1:n_hol; k=k+1;
    pause(0.1); sMovie(k) = getframe(gcf);
end



% x_min = linspace(xlim_zm1(1),xlim_zm2(1),n_zm2);
% x_max = linspace(xlim_zm1(2),xlim_zm2(2),n_zm2);
% 
% y_min = linspace(ylim_zm1(1),ylim_zm2(1),n_zm2);
% y_max = linspace(ylim_zm1(2),ylim_zm2(2),n_zm2);
% 
% 
% for i = 1:n_zm2; k=k+1;
%     set(gca,'xlim',[x_min(i),x_max(i)])
%     set(gca,'ylim',[y_min(i),y_max(i)])
%     pause(0.1); sMovie(k) = getframe(gcf);
% end

for i = 1:n_hol; k=k+1;
    pause(0.1); sMovie(k) = getframe(gcf);
end

iStep_mid = 178;

for i = 1:iStep_mid; k=k+1;
    
    % load ck. coord.
    load([path1,'var_Crack',num2str(i,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a1 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-k','linewidth',1.5);
    end
    
    
    % load ck. coord.
    load([path2,'var_Crack',num2str(i,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a2 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-r','linewidth',1);
    end
    
    
    % load ck. coord.
    load([path1,'var_Crack',num2str(0,'_%05d')]);
    
    for j = 1:length(cCkCrd)
        a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',2);
    end
    
    
%     xlim_zm1 = xlim_zm1 + speed;
    set(gca,'xlim',xlim_zm1,'ylim',ylim_zm1)
    legend([a1(1),a2(1)],'max-hoop','energy min.')
    
    pause(0.01); sMovie(k) = getframe(gcf);
    
    
end

% load ck. coord.
    load([path1,'var_Crack',num2str(628,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a1 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-k','linewidth',1);
    end
    
    
    % load ck. coord.
    load([path2,'var_Crack',num2str(iStep,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a2 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-r','linewidth',1);
    end
    
    
    % load ck. coord.
    load([path1,'var_Crack',num2str(0,'_%05d')]);
    
    for j = 1:length(cCkCrd)
        a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',2);
    end

while xlim_zm1(1) > 0; k=k+1;
    
    
    xlim_zm1 = xlim_zm1 + speed;
    set(gca,'xlim',xlim_zm1,'ylim',ylim_zm1)
    legend([a1(1),a2(1)],'max-hoop','energy min.')
    
    pause(0.01); sMovie(k) = getframe(gcf);
    
    
end

% x_min = linspace(xlim_zm1(1),xlim_nml(1),n_zm1);
% x_max = linspace(xlim_zm1(2),xlim_nml(2),n_zm1);
% 
% y_min = linspace(ylim_zm1(1),ylim_nml(1),n_zm1);
% y_max = linspace(ylim_zm1(2),ylim_nml(2),n_zm1);
% 
% for i = 1:n_zm1; k=k+1;
%     
%     set(gca,'xlim',[x_min(i),x_max(i)])
%     set(gca,'ylim',[y_min(i),y_max(i)])
%     
%     pause(0.1); sMovie(k) = getframe(gcf);
%     
% end

tmp = clock; tmp = ['_',num2str(tmp(4)),'_',num2str(tmp(5))];
movie2avi(sMovie(1:k),[pathS,'mov_CrackGrowth',tmp],...
    'fps',fps,'compression','none')

% clear var.
clear sMovie
% close(h)

%==========================================================================