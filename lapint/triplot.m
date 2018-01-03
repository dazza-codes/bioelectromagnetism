function [hs, hc, contour] = triplot(pnt, dhk, val, mode, levels)

% TRIPLOT make 2D or 3D plot of a triangulated surface and interpolated values
% the surface can be displayed with linear interpolated values or with
% linear interpolated contours of a potential distribution 
% 
% TRIPLOT(pnt, dhk, value, mode, levels) makes a plot of value on the
% surface described by triangles dhk with vertices pnt
%
% dhk can be [], then a delaunay triangulation will be computed (only 2D)
% 
% the visualization mode can be: 
% 	'faces'		plot white triangles only    (value can be [])
%	'face_index'	plot index of each triangle  (value can be [])
% 	'nodes'		plot black vertices only     (value can be [])
%	'node_index'	plot index of each vertex    (value can be [])
% 	'edges'		plot black edges only        (value can be [])
% 	'surface'	make interpolated plot of value on surface
% 	'contour'	make interpolated contour plot of value on surface
% 	'contour_bw'	make interpolated black-white contour plot
% 	'contour_rb'	make interpolated contour plot with red-blue 
% 
% with the optional levels, you can specify the levels at which contours will 
% be plotted 
% 
% See also READ_TRI, READ_BND, ELPLOT, ELPROJ (all in ~roberto/matlab)
% See also PATCH, COLORMAP, VIEW (general Matlab commands)

% updated on Mon Jul 23 12:41:44 MET DST 2001

% (c) 2001 Robert Oostenveld

% start with empty return values

hs      = [];
hc      = [];
contour = [];


figure


% everything is added to the current figure
holdflag = ishold;
hold on

% check the input variables

if isempty(dhk)
  % no triangulation was specified
  if size(pnt,2)==2
    % make a 2d triangulation of the points using delaunay
    dhk = delaunay(pnt(:,1), pnt(:,2));
  elseif strmatch(mode, 'nodes') | strmatch(mode, 'node_index')
    warning('no triangulation specified for 3D vertices')
  else
    error('no triangulation specified for 3D vertices')
  end
end

if nargin<3
  warning('only displaying triangle edges')
  mode='edges';
  val = [];
elseif nargin<4
  warning('default displaying surface')
  mode='surface';
elseif nargin<5
  % determine contour levels
  if ~isempty(val)
    absmax = max(abs([min(val) max(val)]));
    levels = linspace(-absmax,absmax,21);
  else
    levels = [];
  end
end




if strmatch(mode, 'faces')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the faces of the 2D or 3D triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hs = patch('Vertices', pnt, 'Faces', dhk);
  set(hs, 'FaceColor', 'white');
  set(hs, 'EdgeColor', 'none');
end			% mode faces

if strmatch(mode, 'face_index')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the triangle indices (numbers) at each face
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for face_indx=1:size(dhk,1)
    str = sprintf('%d', face_indx);
    tri_x = (pnt(dhk(face_indx,1), 1) +  pnt(dhk(face_indx,2), 1) +  pnt(dhk(face_indx,3), 1))/3;
    tri_y = (pnt(dhk(face_indx,1), 2) +  pnt(dhk(face_indx,2), 2) +  pnt(dhk(face_indx,3), 2))/3;
    tri_z = (pnt(dhk(face_indx,1), 3) +  pnt(dhk(face_indx,2), 3) +  pnt(dhk(face_indx,3), 3))/3;
    h   = text(tri_x, tri_y, tri_z, str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    hs  = [hs; h];
  end 
end

if strmatch(mode, 'edges')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the edges of the 2D or 3D triangulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hs = patch('Vertices', pnt, 'Faces', dhk);
  set(hs, 'FaceColor', 'none');
  set(hs, 'EdgeColor', 'black');
end			% mode edges

if strmatch(mode, 'nodes')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the nodes (vertices) only as points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if size(pnt, 2)==2
    plot(pnt(:,1), pnt(:,2), 'k.');
  else
    plot3(pnt(:,1), pnt(:,2), pnt(:,3), 'k.');
  end
end			% mode edges

if strmatch(mode, 'node_index')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the vertex indices (numbers) at each node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  for node_indx=1:size(pnt,1)
    str = sprintf('%d', node_indx);
    if size(pnt, 2)==2
      h   = text(pnt(node_indx, 1), pnt(node_indx, 2), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    else
      h   = text(pnt(node_indx, 1), pnt(node_indx, 2), pnt(node_indx, 3), str, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    hs  = [hs; h];
  end 
end

if strmatch(mode, 'surface')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot a 2D or 3D triangulated surface with linear interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  hs = patch('Vertices', pnt, 'Faces', dhk, 'FaceVertexCData', val);
  set(hs, 'FaceColor', 'interp');
  set(hs, 'EdgeColor', 'none');
end			% mode surface

if strmatch('contour', mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make a 2D or 3D plot using contours
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

triangle_val = val(dhk);
triangle_min = min(triangle_val, [], 2);
triangle_max = max(triangle_val, [], 2);

for cnt_indx=1:length(levels)
  cnt = levels(cnt_indx);
  use = cnt>=triangle_min & cnt<=triangle_max;
  counter = 0;
  intersect1 = [];
  intersect2 = [];

  for tri_indx=find(use)'
    pos  = pnt(dhk(tri_indx,:), :);
    v(1) = triangle_val(tri_indx,1);
    v(2) = triangle_val(tri_indx,2);
    v(3) = triangle_val(tri_indx,3);
    la(1) = (cnt-v(1)) / (v(2)-v(1));	% abcissa between vertex 1 and 2
    la(2) = (cnt-v(2)) / (v(3)-v(2));	% abcissa between vertex 2 and 3
    la(3) = (cnt-v(3)) / (v(1)-v(3));	% abcissa between vertex 1 and 2
    abc(1,:) = pos(1,:) + la(1) * (pos(2,:) - pos(1,:));
    abc(2,:) = pos(2,:) + la(2) * (pos(3,:) - pos(2,:));
    abc(3,:) = pos(3,:) + la(3) * (pos(1,:) - pos(3,:));
    counter = counter + 1;
    sel     = find(la>=0 & la<=1);
    intersect1(counter, :) = abc(sel(1),:);
    intersect2(counter, :) = abc(sel(2),:);
  end

  % remember the details for external reference
  contour(cnt_indx).level = cnt;
  contour(cnt_indx).n     = counter;
  contour(cnt_indx).intersect1 = intersect1;
  contour(cnt_indx).intersect2 = intersect2;
end

% collect all different contourlevels and plot them
intersect1 = [];
intersect2 = [];
cntlevel   = [];
for cnt_indx=1:length(levels)
  intersect1 = [intersect1; contour(cnt_indx).intersect1];
  intersect2 = [intersect2; contour(cnt_indx).intersect2];
  cntlevel   = [cntlevel; ones(contour(cnt_indx).n,1) * levels(cnt_indx)];
end
   
X = [intersect1(:,1) intersect2(:,1)]';
Y = [intersect1(:,2) intersect2(:,2)]';
C = [cntlevel(:)     cntlevel(:)]';

if size(pnt,2)>2
  Z = [intersect1(:,3) intersect2(:,3)]';
else
  Z = zeros(2, length(cntlevel));
end

hc = [];
for i=1:length(cntlevel)
  
  % make colour contours (default)
  h1 = patch('XData', X(:,i), 'Ydata', Y(:,i), ...
             'ZData', Z(:,i), 'CData', C(:,i), ...
             'facecolor','none','edgecolor','flat', ...
             'linestyle', '-', 'linewidth', 1, ...
             'userdata',cntlevel(i));
  
  if strcmp(mode, 'contour_rb'),
      % make red-blue contours
      if     cntlevel(i)>0,          set(h1,'edgecolor','red');
      elseif cntlevel(i)<0,          set(h1,'edgecolor','blue');
      elseif isequal(cntlevel(i),0), set(h1,'edgecolor','black');
      end
  elseif strcmp(mode, 'contour_bw'),
      % make black-white contours
      set(h1,'edgecolor','black');
      if     cntlevel(i)>0,          set(h1,'edgecolor','black','linestyle', '-');
      elseif cntlevel(i)<0,          set(h1,'edgecolor',[.5 .5 .5],'linestyle', '--');
      elseif isequal(cntlevel(i),0), set(h1,'edgecolor','black','linewidth', 2);
      end
  end
  
  hc = [hc; h1];
end


end			% mode contour

axis vis3d
axis off

if nargout==0
 clear contour hc hs
end

if ~holdflag
  hold off
end

return
