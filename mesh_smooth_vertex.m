function [FV] = mesh_smooth_vertex(FV,nIter),

% mesh_smooth_vertex - smooth a triangulation mesh
%
% [FV] = mesh_smooth_vertex(FV,nIter)
%
% nIter = 10;
%
% See Smith, S. (2002), Fast robust automated brain extraction.  Human Brain
% Mapping, 17: 143-155.  This function uses update component 1 and 2,
% with no contraints for MRI image intensity (update component 3).
%


if ~exist('nIter','var'),
  nIter = 100;
end
if isempty(nIter),
  nIter = 100;
end


[normals,unit_normals] = mesh_vertex_normals(FV);

edges = mesh_edges(FV);


% Define limits for radius of curvature 
% empirically optimized per Smith (2002, fig 6).
rMin =  3.33; % mm
rMax = 10.00; % mm
              % Sigmoid function parameters,
              % "where E and F control the scale and offset of the sigmoid"
E = ((1/rMin) + (1/rMax))/2;
F = 6 /((1/rMin) - (1/rMax));


Nvert = size(FV.vertices,1);

for iteration = 1:nIter,
  
  for index = 1:Nvert,
    
    v = FV.vertices(index,:);
    
    unorm = unit_normals(index,:);
    
    % Find neighbouring vertex coordinates
    vi = find(edges(index,:));
    n = FV.vertices(vi,:);
    nd = edges(index,vi);  % edge distances
    
    % Find neighbour mean location; this is
    % 'mean position of A and B' in
    % figure 4 of Smith (2002)
    nm = mean(n,1);
    
    % Find difference in distance between the vertex
    % of interest and the mean location of its
    % neighbours; this provides 's'
    % in figure 4 of Smith (2002)
    s = nm - v; % toward neighbour mean location
    
    % extract components of s:
    % first the normal component
    sn = dot(s, unorm) * unorm;
    % now the tangent to the local surface
    st = s - sn; % this is an absolute value
    
    % ---------------------------------------
    % movement update components:
    
    % u1 = 1/2 the tangential component
    u1 = st/2;
    
    % u2 is a non-linear function of local surface curvature, where the
    % the local curvature is related to the local radius of curvature,
    % which is given by: r = L^2 / (2 * |sn|);
    % where L is the mean distance of the central vertex to each neighbour
    
    % Find distances between vertex and neighbours, using edge lengths.
    % The mean value is l in Smith (2002, eq. 4)
    L = mean(nd);

    % Calculate radius of local curvature, solve Smith (2002, eq. 4)
    r = (L^2) / (2 * vector_magnitude(sn,v));

    % sigmoid function of r is used for update fraction of sn
    f2 = (1 + tanh(F * (1/r - E))) / 2;
    
    u2 = f2 * sn;
    
    u = u1 + u2;
    
    vu = v + u;
    
    % update the FV data
    FV.vertices(index,:) = vu;
    
  end
end

return
