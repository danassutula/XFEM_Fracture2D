%
% Plot stress intensity factors K1 and K2 
% for a particular crack tip vs. crack length
%
% need to define:
%   plot_SIF_crack
%   plot_SIF_cktip
%

while ~exist('plot_SIF_crack','var') || isempty(plot_SIF_crack) || ...
      ~exist('plot_SIF_cktip','var') || isempty(plot_SIF_cktip)
    
    plot_SIF_crack = input('plot_SIF_crack = ');
    plot_SIF_cktip = input('plot_SIF_cktip = ');
    
end

i_crk = plot_SIF_crack;
i_tip = plot_SIF_cktip;

if SIF_mode1{1}(i_crk,i_tip) ~= 0 
    fprintf(['\nPlotting SIF''s ', ...
        'for crack i_crk = %i ', ...
        'tip i_tip = %i\n'], ...
        i_crk,i_tip)
else
    fprintf(['\nThe SIF''s ', ...
        'for crack i_crk = %i ', ...
        'tip i_tip = %i are undefined\n'], ...
        i_crk,i_tip)
end

scrsz = get(0,'ScreenSize'); % screen position
scrsz = scrsz([3,4])-scrsz([1,2])+1; % screen size

if ~exist('fig_crk','var') || ~ishandle(fig_crk)
    fig_crk = figure;
else
    figure(fig_crk)
end

if exist('szfig_crk','var') 
    set(fig_crk,'Color','w','OuterPosition',szfig_crk)
else
    set(fig_crk,'Color','w','OuterPosition', ...
        [scrsz(1)/2+1,scrsz(2)/2+1,scrsz(1)/2,scrsz(2)/2]) % left-top
end

hold on;
axis equal;

PlotDomain
PlotCracks

plot(cCkCrd_trm{i_crk}(:,1),cCkCrd_trm{i_crk}(:,2),'-b', ...
    'markerfacecolor','w','markeredgecolor','b',...
    'linewidth',1.5,'markersize',5)

switch i_tip
    case 1
        plot(cCkCrd_trm{i_crk}(1,1),cCkCrd_trm{i_crk}(1,2),...
            'ob','markerfacecolor','b')
    otherwise % 2
        plot(cCkCrd_trm{i_crk}(end,1),cCkCrd_trm{i_crk}(end,2),...
            'ob','markerfacecolor','b')
end

if iStep == nStep && mTpRdi(i_crk,i_tip) > 0
    % finished post-processing last solution increment
    q = 1:iStep+1;
else
    % last material state was not post-processed; could be due:
    % (1) a large energy change or (2) negative strain energy
    q = 1:iStep;
end

tpsif = zeros(q(end),2);
cklen = zeros(q(end),1);

for ii = 1:q(end)
    
    tpsif(ii,1) = SIF_mode1{ii}(i_crk,i_tip);
    tpsif(ii,2) = SIF_mode2{ii}(i_crk,i_tip);
    cklen(ii)   = cCkLns{ii}(i_crk);
    
end

if ~exist('fig_sif','var') || ~ishandle(fig_sif)
    fig_sif = figure;
else
    figure(fig_sif)
end

clf(fig_sif); hold on;
set(fig_sif,'Color','w','OuterPosition', ...
    [scrsz(1)/2+1,1,scrsz(1)/2,scrsz(2)/2]) % right-bottom

h_plot = [0,0];
h_plot(1) = plot(cklen,tpsif(:,1),'-or',...
    'linewidth',1,'markerfacecolor','w'); % K1
h_plot(2) = plot(cklen,tpsif(:,2),'-sb',...
    'linewidth',1,'markerfacecolor','w'); % K2
    
if length(cklen) > 1
    set(gca,'xlim',[cklen(1),cklen(end)])
end

legend(h_plot,'K_I','K_{II}')
ylabel('SIF''s: mode-I and mode-II')
xlabel(sprintf('Length of crack (#%i)',i_crk))

title({'Crack tip stress intensity factors vs. crack length',...
    ['(for crack i\_crk = ',num2str(i_crk), ...
    ' and crack tip i\_tip = ',num2str(i_tip),')']})

clear i_crk i_tip
clear tpsif cklen

%--------------------------------------------------------------------------
