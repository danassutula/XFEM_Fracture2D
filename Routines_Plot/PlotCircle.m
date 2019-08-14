function Plot_Circle(r,x0,y0,linespec)

gcf; hold on;

n = 1e3;

b = linspace(0,2*pi,n);
b(end+1) = b(1);

x = x0+r*cos(b);
y = y0+r*sin(b);

if nargin == 3
   linespec = 'k'; 
end
    
plot(x,y,linespec) % ,'linewidth',1.4

end