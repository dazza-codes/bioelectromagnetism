function [potential] = mesh_interp(mesh, potential, depth),

% mesh_interp - Computes a linear interpolation
% 
% Useage: [potential] = mesh_interp(mesh, index, potential, depth)
% 
% This function calculates a linear interpolation of potential values at all
% unknown vertices of a mesh, given known potential values at a subset of
% the mesh vertices.  The function assumes that potential is the same length
% as mesh.vertices, with unknown values at locations given by index.
% 
% mesh.vertices   - vertices of mesh, Nx3 Cartesian XYZ
% mesh.faces      - triangulation of vertices
%
% potential - an Nx3 vector of scalar potentials on the mesh surface.
%
% depth - A parameter to define the depth of vertex connectivity for
% the interpolation (default = 1 for nearest neighbours).
%

% $Revision: 1.1 $ $Date: 2008/05/04 18:25:00 $

% Licence:  GNU GPL, no implied or express warranties
% History:  12/2007, (c) Darren L. Weber
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


error('This function in development.')

if ~exist('depth','var'),
  depth = 1;
end
if isempty(depth),
  depth = 1;
end

% -- Get the mesh connectivity in mesh.connect
mesh = mesh_connect(mesh);

% -- Calculate edge lengths of triangulation into mesh.edge
mesh.edge = mesh_edges(mesh);

% -- Get the index values of the known and unknown potential values

[KnownIndex, i, i] = union(index,index);
if ~isequal(length(KnownIndex),length(index)),
    warning('Trimming duplicate values from index\n');
end
keepindex = sort(i);
repindex = setdiff(1:length(index),sort(i));

KnownIndex = index(keepindex); % unsort KnownIndex
clear index

% find 'unknown' indices
UnknownIndex = setdiff(1:size(mesh.vertices,1), KnownIndex);

for i = 1:length(UnknownIndex),
  
  v = UnknownIndex(i);
  
  neighbourIndex = v;
  for d = 1:depth,
    [x, neighbourIndex] = find(FV.connect(neighbourIndex,:));
    neighbourIndex = unique(neighbourIndex);
    % exclude UnknownIndex values from the interpolation
    neighbourIndex = setdiff(UnknownIndex, neighbourIndex);
  end
  
  % Calculate weights for each neighbour, by distance to v.
  
  
  
  if length(neighbourIndex),
    neighbourValue = potential(neighbourIndex);
    
    % multiply values by weights
    
    potential(v) = mean(neighbourValue);
  end
  
end

return
