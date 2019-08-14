
%==========================================================================
% Plot Cracks (from file)
%==========================================================================
%%

close all

% Figure formating
run([cd,'\other_tools\FigFormat.m'])

path0 = [cd,'\Results\wedgeSplitting_compareEnergyMin_doubleCantilever'];

open([path0,'\domain_3x1.fig']);
hold on; pause(0.1); fps = 8;
FigResize(szfig_720p,szfnt_big)

title('Fracture process');
xlabel('x (mm)');
ylabel('y (mm)');

path_exF = [path0,'\explicit_msh=600x1800\var\'];
path_exp = [path0,'\explicit_msh=300x900\var\'];
path_imp = [path0,'\implicit_msh=300x900\var\'];
path_avg = [path0,'\implicit_msh=300x900_wAvg\var\'];
path_sve = [path0,'\']; % saving path

nStep_exF = 280; % exF
nStep_exp = 136; % exp
nStep_imp = 137; % imp
nStep_avg = 140; % avg

iStep = nStep_exp

k     = 0;

speed = -0.01;

n_zm1 = fps * 1;
n_zm2 = fps * 1;
n_hol = fps * 0.5;

xlim_nml = get(gca,'xlim');
ylim_nml = get(gca,'ylim');

xlim_zm1 = [ 2.342991694678656   3.044866775974183];
ylim_zm1 = [-0.138449807294693   0.138046436852030];

% load ck. coord.
load([path_exp,'var_Crack',num2str(0,'_%05d')]);
% cCkCrd = sCrack.cCkCrd;

for j = 1:length(cCkCrd)
    a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',2);
end
% legend(a0,'pre-existing flaws')


x_min = linspace(xlim_nml(1),xlim_zm1(1),n_zm1);
x_max = linspace(xlim_nml(2),xlim_zm1(2),n_zm1);

y_min = linspace(ylim_nml(1),ylim_zm1(1),n_zm1);
y_max = linspace(ylim_nml(2),ylim_zm1(2),n_zm1);

for i = 1:n_hol; k=k+1;
    pause(0.01); sMovie(k) = getframe(gcf);
end

for i = 1:n_zm1; k=k+1;
    set(gca,'xlim',[x_min(i),x_max(i)])
    set(gca,'ylim',[y_min(i),y_max(i)])
    pause(1e-3); sMovie(k) = getframe(gcf);
end

for i = 1:n_hol; k=k+1;
    pause(1e-3); sMovie(k) = getframe(gcf);
end

% iStep_imp = 137;
for i = 1:100; k=k+1;
    
    % load ck. coord.
    load([path_exF,'var_Crack',num2str(i,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a1f = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'.k','linewidth',1);
    end
    
    
    % load ck. coord.
    load([path_exp,'var_Crack',num2str(0,'_%05d')]);
    
    for j = 1:length(cCkCrd)
        a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',1.1);
    end
    
    
%     xlim_zm1 = xlim_zm1 + speed;
    set(gca,'xlim',xlim_zm1,'ylim',ylim_zm1)
    legend([a1f(1)],'max-hoop ( >2mln. DOF )')
    
    pause(1e-3); sMovie(k) = getframe(gcf);
    
    
end

for i = 1:n_hol; k=k+1;
    pause(1e-3); sMovie(k) = getframe(gcf);
end




% load ck. coord.
load([path_exF,'var_Crack',num2str(nStep_exF,'_%05d')]);
cCkCrd = sCrack.cCkCrd;

for j = 1:length(cCkCrd)
    a1f = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'.k','linewidth',1);
end

% load ck. coord.
load([path_exp,'var_Crack',num2str(nStep_exp,'_%05d')]);
cCkCrd = sCrack.cCkCrd;

for j = 1:length(cCkCrd)
    a1 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-r','linewidth',1.5);
end

% load ck. coord.
load([path_imp,'var_Crack',num2str(nStep_imp,'_%05d')]);
cCkCrd = sCrack.cCkCrd;

for j = 1:length(cCkCrd)
    a2 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-b','linewidth',1.5);
end

for i = 1:n_hol; k=k+1;
    pause(1e-3); sMovie(k) = getframe(gcf);
end




for i = 1:iStep; k=k+1;
    
%     % load ck. coord.
%     load([path_exp,'var_Crack',num2str(i,'_%05d')]);
%     cCkCrd = sCrack.cCkCrd;
%     
%     for j = 1:length(cCkCrd)
%         a1 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-r','linewidth',1.5);
%     end
%     
%     
%     % load ck. coord.
%     load([path_imp,'var_Crack',num2str(i,'_%05d')]);
%     cCkCrd = sCrack.cCkCrd;
%     
%     for j = 1:length(cCkCrd)
%         a2 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-b','linewidth',1.5);
%     end
    
    
    % load ck. coord.
    load([path_avg,'var_Crack',num2str(i,'_%05d')]);
    cCkCrd = sCrack.cCkCrd;
    
    for j = 1:length(cCkCrd)
        a3 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-g','linewidth',1.5);
    end
    
    
    % load ck. coord.
    load([path_exp,'var_Crack',num2str(0,'_%05d')]);
    
    for j = 1:length(cCkCrd)
        a0 = plot(cCkCrd{j}(:,1),cCkCrd{j}(:,2),'-w','linewidth',2);
    end
    
    
    xlim_zm1 = xlim_zm1 + speed;
    set(gca,'xlim',xlim_zm1,'ylim',ylim_zm1)
    leg = legend([a1f(1),a1(1),a2(1),a3(1)],'max-hoop    (>2.0mln. DOF)',...
                                            'max-hoop    (\approx0.5mln. DOF)',...
                                            'energy min. (\approx0.5mln. DOF)',...        
                                            'improved     (\approx0.5mln. DOF)')
    set(leg,'fontsize',16);
                                        
    pause(1e-3); sMovie(k) = getframe(gcf);
    
    
end

for i = 1:n_hol; k=k+1;
    pause(1e-3); sMovie(k) = getframe(gcf);
end

x_min = linspace(xlim_zm1(1),xlim_nml(1),n_zm1);
x_max = linspace(xlim_zm1(2),xlim_nml(2),n_zm1);

y_min = linspace(ylim_zm1(1),ylim_nml(1),n_zm1);
y_max = linspace(ylim_zm1(2),ylim_nml(2),n_zm1);

for i = 1:n_zm1; k=k+1;
    
    set(gca,'xlim',[x_min(i),x_max(i)])
    set(gca,'ylim',[y_min(i),y_max(i)])
    
    pause(1e-3); sMovie(k) = getframe(gcf);
    
end

tmp = clock; tmp = ['_',num2str(tmp(4)),'_',num2str(tmp(5))];
movie2avi(sMovie(1:k),[path_sve,'mov_CrackGrowth',tmp],...
    'fps',fps,'compression','none')

% clear var.
clear sMovie
% close(h)

%==========================================================================