function d = ElemQuality_t3(x)
%--------------------------------------------------------------------------
%
% GeoElm_quality:
%   Trinagle quality measure based on an analysis of interpolation error
%   (Bank and Smith): the area devided by the sum of squared edge lengths.
%   P.S.: this measure slighlty favors sharp triangles over flat triangles
%
% INPUT:
%   x (3x2) - coordinates of element (triangle) vertices
% OUTPUT:
%   d (1x1) - a measure of the distortion of a trinagle, d \in [0 .. 1]
%
%--------------------------------------------------------------------------

% relative coordinates of element vertices
x = [x(2,:)-x(1,:);x(3,:)-x(1,:);x(3,:)-x(2,:)];

% triangualr element quality measure (Bank and Smith)
d = 0.5*(x(1,1)*x(2,2)-x(1,2)*x(2,1)) / sum(x(:,1).^2+x(:,2).^2);

% normalize
d = 6.928203230275510*d;

end