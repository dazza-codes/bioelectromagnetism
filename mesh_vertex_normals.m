function [vertNormals,vertNormalsUnit] = mesh_vertex_normals(FV,index,area)

% mesh_vertex_normals - Calculate vertex surface normals
% 
% [vertNormals,vertNormalsUnit] = mesh_vertex_normals(FV,index,area)
% 
% FV.vertices - vertices of mesh, Nx3 Cartesian XYZ
% FV.faces    - triangulation of vertices (Mx3 matrix)
% 
% index is an array of vertex indices (default is all vertices)
% 
% area = 0, vertex normal is the mean of all face normals (default)
% area = 1, vertex normal is the mean of all face normals, weighted by the
% area of each face (the weighting is faceNormal * faceArea).
% 
% normals       - vertex surface normals (Nx3 matrix)
% unit_normals  - normalized normals!
% 
% This routine first calculates the surface normals of each face.  It then
% finds all the faces that contain a specific vertex and takes the mean of
% those face normals to estimate the vertex normal (it can be weighted by
% the face area).
% 
% If the faces are defined according to the right-hand rule, all their
% normals will be "outward" normals and the average should be a sensible
% value for the vertex normal.  See mesh_face_normals for details.
% 
% See also the VertexNormals property of the patch and surface commands.
% 

% $Revision: 1.6 $ $Date: 2007/04/13 19:22:01 $

% 
% Licence:  GNU GPL, no implied or express warranties
% History:  04/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% first get some face data
[faceNormals,faceNormalsUnit,faceCentroids,faceArea] = mesh_face_normals(FV);

tic;
fprintf('...calculating vertex surface normals...');

if ~exist('area','var'),
  area = 0;
end
if isempty(area),
  area = 0;
end

if area,
  % Weight the faceNormals by the faceArea
  faceNormals = faceNormals .* repmat(faceArea,1,3);
end


Nvert = size(FV.vertices,1);

if ~exist('index','var'),
  index = 1:Nvert;
end
if isempty(index),
  index = 1:Nvert;
end

vertNormals = zeros(length(index),3);
vertNormalsUnit = vertNormals;

for i = 1:length(index),
    
    vi = index(i);
    v = FV.vertices(vi,:);
    
    % To calculate a vertex normal you need to perform an iteration beginning
    % with a triangle that holds the vertex and perform a cross product on
    % the two vectors of that triangle that extend from that vertex. This
    % will result in a normal for that triangle. Then go around to each
    % triangle containing that vertex and perform a cross product for that
    % triangle. You will end up with a set of normals for all the triangles
    % (faces). To compute the vertex normal, take the mean of all the face
    % normal vectors.
    
    % get all the faces that contain this vertex
    [faceIndexI,faceIndexJ] = find(FV.faces == vi);
    
    Nface = length(faceIndexI);
    
    faceNorms = faceNormals(faceIndexI,:);
    
    vertNormals(i,:) = sum(faceNorms,1) / Nface;
    
    vnorm = vector_magnitude(vertNormals(i,:));
    vertNormalsUnit(i,:) = vertNormals(i,:) / vnorm;
    
end

t=toc;
fprintf('done (%5.2f sec).\n',t);

return
