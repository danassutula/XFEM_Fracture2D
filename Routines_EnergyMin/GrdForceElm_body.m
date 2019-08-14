function [grd_Fe,gd2_Fe] = GrdForceElm_body(type,dx,x,N,dNdX,W,f)

%--------------------------------------------------------------------------

% Fe = int N*f dOmg

% (!!!) if the element is enriched, the enrichment must remain constant

%--------------------------------------------------------------------------

switch type % derivative type
    case 'inc' % 2nd deriv. are zero
        d2x = zeros(size(dx));
    case 'rot' % 2nd deriv. are non zero
        d2x = [-dx(:,2),dx(:,1)];
    otherwise
        error('Unknown type of differentiation; type = ''inc'' or ''rot''')
end

nJ = size(x,1);
nN = size(N,2);
nD = 2*nN;

grd_Fe = zeros(nD,1);
gd2_Fe = zeros(nD,1);

N1 = zeros(1,nN);
N2 = zeros(1,nN);

% need Jackobian only
dNdX = dNdX(:,1:nJ);

q = 1:2; % id deriv.
for i = 1:length(W)
    
    grd_DJ = dNdX(q(1),:)'*dNdX(q(2),:);
    grd_DJ = grd_DJ - grd_DJ';
    
    del_DJ =  dx(:,1)'*grd_DJ*x(:,2) + x(:,1)'*grd_DJ*dx(:,2);
    dl2_DJ = d2x(:,1)'*grd_DJ*x(:,2) + 2*dx(:,1)'*grd_DJ*dx(:,2) + x(:,1)'*grd_DJ*d2x(:,2);
    
    N1 = N1 + N(i,:)*(del_DJ*W(i));
    N2 = N2 + N(i,:)*(dl2_DJ*W(i));
    
    q = q+2;
    
end

grd_Fe(1:2:end-1) = N1*f(1);
grd_Fe(2:2:end)   = N1*f(2);

gd2_Fe(1:2:end-1) = N2*f(1);
gd2_Fe(2:2:end)   = N2*f(2);

end