
nCrack = 1;
cCkCrd = cell(nCrack,1);

% a_crk = 5.0;
% b_crk = 1.5;

% a_crk = 6.0;
% b_crk = 1.0;

a_crk = 6.0;
b_crk = 2.5;

mCkCrd = zeros(2,2);
mCkCrd(1,:) = [10-a_crk,-1];
mCkCrd(2,:) = [10-a_crk, b_crk];

cCkCrd{1} = mCkCrd;
