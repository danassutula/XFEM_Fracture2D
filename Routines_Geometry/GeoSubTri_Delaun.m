function [mSbCrd,mSbNod] = GeoSubTri_Delaun(mSbCrd,mEdges)

% wrapper function for Delaunay triangulation

if nargin == 1
    mEdges = [];
elseif ~isempty(mEdges) && size(mEdges,2)~=2
    error('mEdges must have size: [n_edges,2]')
end

try % exist('delaunayTriangulation','builtin')
    
    if isempty(mEdges)
        tri = delaunayTriangulation(mSbCrd);
    else
        tri = delaunayTriangulation(mSbCrd,mEdges);
    end
        
    mSbNod = tri.ConnectivityList;
    mSbCrd = tri.Points;
    
catch % exist('DelaunayTri','builtin')
    
    if isempty(mEdges)
        tri = DelaunayTri(mSbCrd);
    else
        tri = DelaunayTri(mSbCrd,mEdges);
    end
    
    mSbNod = tri.Triangulation;
    mSbCrd = tri.X;
    
end

n = size(mSbNod,1); % no. tri.
A = zeros(n,1); % initialise area vector

for i = 1:n
    
    x = mSbCrd(mSbNod(i,:),:); % vertice coord.
    A(i) = (x(2)-x(1))*(x(6)-x(4)) - (x(3)-x(1))*(x(5)-x(4)); % get 2 x area
    
end

mSbNod = mSbNod(A > max(A)*1e-12,:); % keep el. if area above tol.

end