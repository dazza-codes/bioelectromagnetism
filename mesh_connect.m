function FV = mesh_connect(FV)

% mesh_connect - Calculate triangulation connectivity
% 
% FV = mesh_connect(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices
% 
% FV.connect - a sparse matrix of vertex connections by index value.
%
% To find all the vertices connected to any one vertex (v), use:
%
% x = find(FV.connect(v,:));
%
% To obtain their coordinates, use FV.vertices(x,:).  For example:
%
% patch('vertices', FV.vertices, 'faces', FV.faces, 'facecolor', 'r')
% hold on
% plot3(FV.vertices(v,1), FV.vertices(v,2), FV.vertices(v,3), 'bo')
% plot3(FV.vertices(x,1), FV.vertices(x,2), FV.vertices(x,3), 'go')
% rotate3d
%

% $Revision: 1.2 $ $Date: 2007/12/15 01:01:42 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  12/2007, (c) Darren L. Weber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...searching for vertex connectivity...');

nVertex = size(FV.vertices,1);

FV.connect = sparse(nVertex,nVertex);

for v = 1:nVertex,
  
  [i, j] = find(FV.faces == v);
  vert = unique(FV.faces(i,:));
  vert = setdiff(vert, v);
  
  FV.connect(v, vert) = 1;
  
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return
