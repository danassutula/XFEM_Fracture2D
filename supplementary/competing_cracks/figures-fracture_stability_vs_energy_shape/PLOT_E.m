
% stable stable
% g = -[1.5;2.5];
% H = [1,-1.;-1.,1];

% stable unstable
% g = -[2.5;2];
% H = [1,2;2,1];

% unstable stable
% g = -[1.5;0.5];
% H = [-2,-5.;-5.,-1];

% unstable stable
% g = -[1.5;0.5];
% H = [-2,-5.;-5.,-1];

% close all
% clear all

figure; hold on;
set(gca,'color','w')
set(gcf,'color','w')

g = -[1.5;0.5];
H = [-2,-5.;-5.,-1];

C = eye(2)-ones(2)/2;
u_star = eig(C*H*C)


[X,Y] = meshgrid(linspace(0,1,10));
Z = 1/2*diag([X(:),Y(:)]*H*[X(:),Y(:)]');
Z = Z + [X(:),Y(:)]*g; % add gradien of energy

Z = reshape(Z,size(X));
h_psi = mesh(X,Y,Z,'facealpha',0);


x_sol = diag(X);
y_sol = 1-diag(Y);
z_sol = 1/2*diag([x_sol,y_sol]*H*[x_sol,y_sol]');
z_sol = z_sol + [x_sol,y_sol]*g

h_sol = plot3(diag(X),1-diag(Y),z_sol+1.0e-2, ...
    '-ob','linewidth',1,'markerfacecolor','w');

%  {\scriptstyle\sum}
% h_leg = legend([h_sol],'$$\Psi(\dot{a})$$, s.t. $$\sum \dot a = 1$$');
h_leg = legend([h_sol],'$$\mathcal{E}(\Delta \ell_1^k,\Delta \ell_2^k)$$ s.t. $$\Delta \ell^k_1+\Delta \ell^k_2 = \Delta a^k$$',...
    'Location','NorthEast'); % 'best'
set(h_leg,'Interpreter','latex');

hx = xlabel('$$\Delta \ell_1$$','Interpreter','latex');
hy = ylabel('$$\Delta \ell_2$$','Interpreter','latex');
hz = zlabel('$$\mathcal{E}(\Delta \ell_1,\Delta \ell_2)$$','Interpreter','latex');

% title('$$\Psi(\dot{a}) = {\frac{1}{2}} {H_s}_{ij} a_i a_j$$','Interpreter','latex');


view(115,20)
grid on


scrsz = get(0,'ScreenSize'); scrsz = scrsz([3,4]);
szfig = [scrsz(1)/2-280,scrsz(2)/2-210,560,420]; % std. fig. size

% szfnt_std = 11; FigResize(szfig,szfnt_std,gcf);
% szfnt_big = 22; FigResize(szfig,szfnt_big,gcf);
szfnt_big = 22; FigResize(szfig,szfnt_big,gcf);
set(h_leg,'fontsize',16)

% shading faceted
% shading interp


% SavePDF(['plot(E)_result(Hs-Vs-Hs_star)'],gcf)
% SavePDF(['plot(E)_grwoth(unstable)_front(stable)'],gcf)

