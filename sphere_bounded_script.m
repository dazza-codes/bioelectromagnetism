
sphere = sphere_tri('ico',4,4);

Xbound = 2;
Ybound = 2;
Zbound = 2;

for i = 1:size(sphere.vertices,1),
  
  x = sphere.vertices(i,1);
  y = sphere.vertices(i,2);
  z = sphere.vertices(i,3);
  
  if x > Xbound, delete_index(end) = i; end
  if y > Ybound, delete_index(end) = i; end
  if z > Zbound, delete_index(end) = i; end
  
  % at this point you know if it lies outside
  % a bounding box, so try to find the distance
  % of this point from the bounding box and
  % determine whether it is less than the
  % length of a triangle edge.  If it is, then
  % maybe don't delete this vertex, just place 
  % it on the edge of the box.  This is the
  % tricky bit, as what you might try to do
  % is maintain the integrity of the
  % triangulation at the bounding box edges.
  
  vertex2boxDistance = 0; % you need to calc this
  vertexEdgeDistance = 0; % you need to calc this
  
  if vertex2boxDistance < vertexEdgeDistance,
    % don't delete the vertex, just move it
    % consider using dsearchn function for this
  end
end

