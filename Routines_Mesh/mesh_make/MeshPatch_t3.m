function [p,t] = MeshPatch_t3(fd,h0,ptfix)

% INPUT:
%   fd      - level set function handle
%   h0      - desired element size
%   ptfix   - boundary points (closed) given in a consecutive order
% OUTPUT:
%   p       - mesh node coordinates
%   t       - element connectivities

% OUTLINE:
% start by creating a grid of equilateral triangales; discard points
% outside poly-line boundary and points close to the boundary (<h0);
% perform delaunay triangulation; discard zero-area elements; discard
% elements lying outside poly-line boundary; fix element node order.

if ptfix(1,:) ~= ptfix(end,:)
    error('End nodes must correspond for a closed polyline boundary')
end

% center of the polygon
pc = sum(ptfix,1)/size(ptfix,1);

% vertical spacing
hv = h0*sqrt(3)/2;

% make initial distribution within polygon (equilateral tri.)
[x,y] = meshgrid([sort(pc(1)-h0:-h0:min(ptfix(:,1))),pc(1):h0:max(ptfix(:,1))],...
    [sort(pc(2)-hv:-hv:min(ptfix(:,2))),pc(2):hv:max(ptfix(:,2))]);

x(2:2:end,:) = x(2:2:end,:)+h0/2;
p = [x(:),y(:)]; % node coord.

% rmv. pt. outside
[H,S] = fd(p,ptfix);
p = p(H>0 & S>h0*0.5,:); % (no obtuse triangles, e.g. >90deg.)

% get unq. ptfix and append new pt.
p = [ptfix(1:end-1,:);p];

% triangulate
t = delaunayn(p);

% to discard zero-area el.
A = zeros(size(t,1),1);

for i = 1:length(A); x = p(t(i,:),:); % get 2x area
    A(i) = (x(2)-x(1))*(x(6)-x(4))-(x(3)-x(1))*(x(5)-x(4));
end

% keep el. if area above tol.
t = t(A > h0^2*1e-9,:);

% to discard el. outside ptfix
x = zeros(size(t,1),2);

for i = 1:size(x,1)
    x(i,:) = sum(p(t(i,:),:),1)/3;
end

% rmv. el. outside
H = fd(x,ptfix);
t = t(H>0,:);

% anti-clockwise nd. order
t = GeoElm_ndord(t,p);

end