function Ke = MtxStiffElm(x_elm,D,dNdX,W)

nn = size(x_elm,1);
nd = 2*size(dNdX,2);

Ke = zeros(nd);
B  = zeros(3,nd);
iJ = zeros(2);

qJ = 1:nn;
qy = 2:2:nd;
qx = qy-1;
qd = 1:2;

for w = W(:)'
    
    J  = dNdX(qd,qJ)*x_elm;
    DJ = J(1)*J(4)-J(2)*J(3);
    
    iJ(1) =  J(4)/DJ; iJ(3) = -J(3)/DJ;
    iJ(2) = -J(2)/DJ; iJ(4) =  J(1)/DJ;
    
    dNdx = iJ*dNdX(qd,:); qd = qd+2;
    
    B(1,qx) = dNdx(1,:);
    B(3,qy) = dNdx(1,:);
    B(2,qy) = dNdx(2,:);
    B(3,qx) = dNdx(2,:);
    
    Ke = Ke + B'*(D*(DJ*w))*B;
    
end
end
