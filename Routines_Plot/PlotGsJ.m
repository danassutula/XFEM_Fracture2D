%
% Plot the J-integral for a particular crack tip vs. crack length
%
% need to define:
%   plot_crack
%   plot_cktip
%

while ~exist('plot_GsJ_crack','var') || isempty(plot_GsJ_crack) || ...
      ~exist('plot_GsJ_cktip','var') || isempty(plot_GsJ_cktip)
    
    plot_GsJ_crack = input('plot_GsJ_crack = ');
    plot_GsJ_cktip = input('plot_GsJ_cktip = ');
    
end

i_crk = plot_GsJ_crack;
i_tip = plot_GsJ_cktip;

if JIntegral{1}(i_crk,i_tip) ~= 0 
    fprintf(['\nPlotting J-integral value ', ...
        'for crack i_crk = %i ', ...
        'tip i_tip = %i\n'], ...
        i_crk,i_tip)
else
    fprintf(['\nThe J-integral value', ...
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

tpGsJ = zeros(q(end),1);
cklen = zeros(q(end),1);

for ii = 1:q(end)
    
    tpGsJ(ii,1) = JIntegral{ii}(i_crk,i_tip);
    cklen(ii)   = cCkLns{ii}(i_crk);
    
end

if ~exist('fig_GsJ','var') || ~ishandle(fig_GsJ)
    fig_GsJ = figure;
else
    figure(fig_GsJ)
end

clf(fig_GsJ); hold on;
set(fig_GsJ,'Color','w','OuterPosition', ...
    [scrsz(1)/2,1,scrsz(1)/2,scrsz(2)/2]) % left-bottom

h_plot = plot(cklen,tpGsJ,'-^m',...
    'linewidth',1,'markerfacecolor','w'); % J
    
if length(cklen) > 1
    set(gca,'xlim',[cklen(1),cklen(end)])
end

legend(h_plot,'J')
ylabel('J-integral value, J')
xlabel(sprintf('Length of crack (#%i)',i_crk))

title({'Crack tip J-integral value vs. crack length',...
    ['(for crack i\_crk = ',num2str(i_crk), ...
    ' and crack tip i\_tip = ',num2str(i_tip),')']})

clear i_crk i_tip
clear tpGsJ cklen

%--------------------------------------------------------------------------
