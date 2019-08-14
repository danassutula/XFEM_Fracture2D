function h = ElemSize(x)

n = size(x,1);
q = [1:n;2:n,1];
s = zeros(n,1);

for i = 1:n
    s(i) = ( x(q(2,i),1)-x(q(1,i),1) )^2 + ( x(q(2,i),2)-x(q(1,i),2) )^2 ;
end

h = sqrt(max(s));
