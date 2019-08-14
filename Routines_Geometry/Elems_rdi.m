function idelm = Elems_rdi(ndcrd,elnod,x_tip,r_tip)

% Elems_tip0:
%   returns elements that have at least one node within the disk r_tip that
%   is centered about x_tip; ndcrd are the nodal coordinates; elnod is the
%   element topology.

idelm = find(any(ismember(elnod,find((ndcrd(:,1)-x_tip(1)).^2+(ndcrd(:,2)-x_tip(2)).^2<r_tip^2)),2));

end