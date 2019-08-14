function yi = InterpLin(xi,x,y)
%==========================================================================

yi = zeros(size(xi));

if ~issorted(x);
    [x,j] = sort(x); y = y(j);
end

j1 = xi <= x(1);
j3 = xi >= x(end);
j2 = ~(j1 | j3);

m = (y(2)-y(1))/(x(2)-x(1));
for i = find(j1(:))'
    yi(i) = y(1) + (xi(i)-x(1))*m;
end

for i = find(j2(:))'
    k1 = find( x < xi(i),1,'last'); k2 = k1 + 1;
    yi(i) = y(k1) + (xi(i)-x(k1))*(y(k2)-y(k1))/(x(k2)-x(k1));
end

m = (y(end)-y(end-1))/(x(end)-x(end-1));
for i = find(j3(:))'
    yi(i) = y(end) + (xi(i)-x(end))*m;
end

%==========================================================================

end