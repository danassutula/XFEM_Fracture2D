
close all
clear all

%%% CHOICES:
% TYPE_GS = 'uniform'
TYPE_GS = 'nonunif'

%%% CHOICES:
TYPE_HS = 'indefn'  % is indefinite
% TYPE_HS = 'posdef' % is positive def.
% TYPE_HS = 'negdef' % is negative def.

%%% CHOICES:
N_TIPS = [3,6,12]

line_width      = 1;
line_style      = {':^r','-.sb','--dk'};
marker_color    = 'w';
marker_size     = 8;

h_plot = zeros(1,length(N_TIPS));
for i_tips = 1:length(N_TIPS)

NAME = sprintf('GsType(%s)_HsType(%s)_nTips(%i)',...
    TYPE_GS,TYPE_HS,N_TIPS(i_tips));

load(['cases-1/results/Gradient-SingleTrial-FixedLength/',NAME,'_results']);

n_mesh = INPUTS.n_mesh;
da_inc = zeros(n_mesh,1);
da_err = zeros(n_mesh,1);

for i = 1:n_mesh
    
    da_h = RESULTS(i).da_h;
    da_e = RESULTS(n_mesh).da_e;
    
    da_err(i) = norm(da_h-da_e)/norm(da_e);
    da_inc(i) = RESULTS(i).da_inc;
    
end

% rescale all lengths to 1
da_inc = da_inc/max(da_inc);

h_plot(i_tips)=loglog(da_inc,da_err,line_style{i_tips},...
    'linewidth',line_width,'markerfacecolor',marker_color,'markersize',marker_size);

hold on;

end

axis tight
xlim = get(gca,'xlim'); xlim(2) = 1;
set(gca,'xlim',xlim);
ylim = get(gca,'ylim'); ylim(:) = [1e-4,10];
set(gca,'ylim',ylim)

h_leg = legend(h_plot,'$n_\mathrm{tips} = 3$',...
    '$n_\mathrm{tips} = 6$','$n_\mathrm{tips} = 12$',...
    'location','southeast'); set(h_leg,'interpreter','latex');

h_x = xlabel('Fracture increment length, $\Delta \ell_\mathrm{inc}/\Delta a_\mathrm{tot}$','interpreter','latex');
% h_y = ylabel('Difference in final crack front, $\sqrt{\textstyle{\sum}_{i=1}^{n_\mathrm{tip}} (\ell^\mathrm{fix}_i - \ell_i)^2 / \textstyle{\sum} \ell_i^2}$');
h_y = ylabel('Difference in final crack front, $\sqrt{\textstyle{\sum}_{i=1}^{n_\mathrm{tip}} \Delta \ell_i^2 / \textstyle{\sum} \ell_i^2}$');

set(gcf,'color','w')
% set(h_y,'unit','character')
set(h_y,'Interpreter','latex');

szfig = [560,420];
szfnt_big = 17; 

FigResize(szfig,szfnt_big,gcf);

NAME = sprintf('GsType(%s)_HsType(%s)_method(%s)',...
    TYPE_GS,TYPE_HS,'GradientSingleTrialFixedLength');

SavePDF(['figures/cases-1/',NAME],gcf)
