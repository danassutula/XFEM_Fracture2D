
% close all
% clear all

figure; hold on;
set(gca,'color','w')
set(gcf,'color','w')


H = [1,3;3,1];
C = eye(2)-ones(2)/2;
u_star = eig(C*H*C)


[X,Y] = meshgrid(linspace(0,1,10));
Z = 1/2*diag([X(:),Y(:)]*H*[X(:),Y(:)]');
Z = reshape(Z,size(X));

h_psi = mesh(X,Y,Z);


x_sol = diag(X);
y_sol = 1-diag(Y);
z_sol = 1/2*diag([x_sol,y_sol]*H*[x_sol,y_sol]');

h_sol = plot3(diag(X),1-diag(Y),z_sol+1.9e-2, ...
    '-ob','linewidth',1,'markerfacecolor','w');


% h_leg = legend([h_sol],'$$\Psi(\dot{a})$$, s.t. $$\sum \dot a = 1$$');
h_leg = legend([h_sol],'$$\Psi(\mathbf{v})$$, s.t. $${\scriptstyle\sum} \mathbf{v} = 1$$');
set(h_leg,'Interpreter','latex');

hx = xlabel('$$v_1$$','Interpreter','latex');
hy = ylabel('$$v_2$$','Interpreter','latex');
hz = zlabel('$$\Psi(\mathbf{v})$$','Interpreter','latex');

% title('$$\Psi(\dot{a}) = {\frac{1}{2}} {H_s}_{ij} a_i a_j$$','Interpreter','latex');


view(-60,15)
grid on


scrsz = get(0,'ScreenSize'); scrsz = scrsz([3,4]);
szfig = [scrsz(1)/2-280,scrsz(2)/2-210,560,420]; % std. fig. size

% szfnt_std = 11; FigResize(szfig,szfnt_std,gcf);
% szfnt_big = 22; FigResize(szfig,szfnt_big,gcf);
szfnt_big = 26; FigResize(szfig,szfnt_big,gcf);

% set(hy,'position',get(hy,'position')+[0,-0.1,0.05]);

% % shading faceted
% shading interp


% SavePDF(['plot(Psi)_result(Hs-Vs-Hs_star)'],gcf)
% SavePDF(['plot(Psi-convex)_result(Hs-concave)'],gcf)
