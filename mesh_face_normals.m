function [faceNormals,faceNormalsUnit,faceCentroids,faceArea] = mesh_face_normals(FV)

% mesh_face_normals - Calculate face surface normals
% 
% [faceNormals,faceNormalsUnit,faceCentroids,faceArea] = mesh_face_normals(FV)
% 
% FV.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces      - triangulation of vertices (Mx3 matrix)
% 
% faceNormals      - face surface normals (Mx3 matrix)
% faceNormalsUnit  - normalized normals!
%
% faceCentroids     - the point at the center of each face
% faceArea          - the face area
%
% When calculating the surface normal, this function assumes the convention
% that the vertices are given in FV.faces in clockwise order:
%
%        1
%       /\
%      /  \
%   e2/    \e1
%    /      \
%   /        \
%  /__________\
% 3	    	   2
%
% We then define edge1 (e1) from vertex1 to vertex2 and edge2 (e2) from
% vertex1 to vertex3.  The cross product of these two edge vectors gives
% the surface normal.  The direction of the normal is either into the
% page or out of the page, depending on the order of the cross product,
%
% edge1 x edge2 = -1 * ( edge2 x edge1 )
%
% So, the order of the vertex points, the direction of the edge vectors
% and their cross products are very important if you want a particular
% direction.
%
% In this function, we assume that the vertices of each face are given in
% clockwise order, when viewed from the "outside".  The resulting surface
% normals are oriented "outward".  Here is an example:
%
% v1 = [ 0 0 0 ]; % origin
% v2 = [ 0 1 0 ]; % y unit vector
% v3 = [ 1 0 0 ]; % x unit vector
% FV.vertices = [v1;v2;v3];
% FV.faces = [ 1 2 3 ];
% [faceNormals,faceNormalsUnit,faceCentroids,faceArea] = mesh_face_normals(FV)
% figure
% Hp = patch('faces',FV.faces,'vertices',FV.vertices,...
%     'facecolor',[.7 .7 .7],'facealpha',0.8,'edgecolor',[.8 .8 .8]); 
% camlight('headlight','infinite'); daspect([1 1 1]); axis vis3d; axis off
% material dull; hold on
% quiver3(faceCentroids(:,1),faceCentroids(:,2),faceCentroids(:,3),...
%     faceNormals(:,1),faceNormals(:,2),faceNormals(:,3));
%


% $Revision: 1.2 $ $Date: 2007/04/13 19:29:31 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;
fprintf('...calculating face surface normals...');

nface = size(FV.faces,1);

faceNormals = zeros(nface,3);
faceNormalsUnit = faceNormals;
faceCentroids = faceNormals;

for f = 1:nface,
    
    vertex_index1 = FV.faces(f,1);
    vertex_index2 = FV.faces(f,2);
    vertex_index3 = FV.faces(f,3);
    
    vertex1 = FV.vertices(vertex_index1,:);
    vertex2 = FV.vertices(vertex_index2,:);
    vertex3 = FV.vertices(vertex_index3,:);
    
    % If the vertices are given in clockwise order, when viewed from the
    % outside, then following calculates the "outward" surface normals.
    
    edge_vector1 = vertex2 - vertex1;
    edge_vector2 = vertex3 - vertex1;
    
    faceNormals(f,:) = cross( edge_vector2, edge_vector1 );
    
    faceNormalsUnit(f,:) = vector_unit(faceNormals(f,:));
    
    % Now find the midpoint between all vertices
    a = (vertex1 + vertex2) ./ 2;
    b = (vertex2 + vertex3) ./ 2;
    c = (vertex3 + vertex1) ./ 2;
    
    % Now find the centroid length of the medians
    faceCentroids(f,:) = b + ( (vertex1 - b) ./3 );
    
end
    
% Area of Triangle = ||AB x AC|| / 2 ; ie, the absolute value of
% the length of the cross product of AB and AC divided by 2
faceNormalsMagnitude = vector_magnitude(faceNormals);
faceArea = abs(faceNormalsMagnitude) ./ 2;

t=toc;
fprintf('done (%5.2f sec).\n',t);

return
