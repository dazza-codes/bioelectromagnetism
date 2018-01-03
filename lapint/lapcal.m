function [lap, edge] = mesh_lapcal(vertex, face)

% MESH_LAPCAL: Compute the laplacian matrix for a triangulation
% 
% Useage:   lap = lapcal(vertices, faces)
% 
% Returns a full matrix that is the Laplacian (2nd spatial
% derivative) of an irregular triangular mesh.  The routine
% is given in:
% 
% Oostendorp T, Oosterom A & Huiskamp G (1989),
% Interpolation on a triangulated 3D surface.
% Journal of Computational Physics, 80: 331-343.
% 
% see also LAPINT
% 

% Licence:  GNU GPL, no implied or express warranties
% History:  (c) 04/2002 Robert Oostenveld/Darren.Weber@flinders.edu.au
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvertex = size(vertex,1);
nface = size(face,1);

% the matrix 'edge' is the connectivity of all vertices
%edge = spalloc(nvertex, nvertex, 3*nface);
edge = zeros(nvertex); % this runs faster
for i=1:nface,
    
    % compute the length of all triangle edges
    edge(face(i,1), face(i,2)) = norm(vertex(face(i,1),:) - vertex(face(i,2),:));
    edge(face(i,2), face(i,3)) = norm(vertex(face(i,2),:) - vertex(face(i,3),:));
    edge(face(i,3), face(i,1)) = norm(vertex(face(i,3),:) - vertex(face(i,1),:));
    % make sure that all edges are symmetric
    edge(face(i,2), face(i,1)) = edge(face(i,1), face(i,2));
    edge(face(i,3), face(i,2)) = edge(face(i,2), face(i,3));
    edge(face(i,1), face(i,3)) = edge(face(i,3), face(i,1));
end

lap = zeros(nvertex);
for i=1:nvertex,
    
    k = find(edge(i,:));        % the indices of the neighbours
    
    ni = length(k);             % the number of neighbours
    
    hi = mean(edge(i,k));       % the average distance to the neighbours
    invhi = mean(1./edge(i,k)); % the average inverse distance to the neighbours
    
    lap(i,i) = -(4/hi) * invhi;
    
    lap(i,k) =  (4/(hi*ni)) * 1./edge(i,k);
    
end

return
