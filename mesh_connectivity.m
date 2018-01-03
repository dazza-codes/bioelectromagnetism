function FV = mesh_connectivity(FV)

% mesh_connectivity - Calculate triangulation connectivity
% 
% FV = mesh_connectivity(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices
% 
% FV.connect    - vertex connections; a sparse matrix
%

% $Revision: 1.1 $ $Date: 2008/05/04 18:25:00 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  07/2002, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...searching for vertex connectivity...');

nVertex = size(FV.vertices,1);

% the 'edge' matrix is the connectivity of all vertices
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
