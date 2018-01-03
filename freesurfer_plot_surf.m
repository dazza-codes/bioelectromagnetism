function [Hf,Hp] = freesurfer_plot_surf(Hf,surf,normals,orientOnly)

% freesurfer_plot_surf
%
% [Hf,Hp] = freesurfer_plot_surf(Hf,surf,normals,orientOnly)
%
% Hf is a figure handle to plot into, using hold on; if empty, a new figure
% is created.
%
% surf is a struct, returned by freesurfer_read_surf, with fields:
%
%    surf.faces     - Tx3 triangles of the surface
%    surf.vertices  - Vx3 vertex coordinates
%
% normals is 0 or 1; if 1, plot the surface normal vectors; see
% freesurfer_read_surf for notes about the order of vertices in the face
% matrix to obtain outward normals. For outward normals, prior to input to
% this function, use:
%
%    surf.faces = surf.faces(:,[1 3 2]);
%
% orientOnly is 0 or 1; if 1, and normals = 1, the normals are converted to
% unit normals, so they only indicate orientation.  This vector field has a
% more uniform appearance.
%
% Hf is a handle to the figure
% Hp is a handle to the patch object
%

% $Revision: 1.4 $ $Date: 2005/08/14 21:22:42 $

% Copyright (C) 2000  Darren L. Weber
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  04/2005, Darren.Weber_at_radiology.ucsf.edu
%           06/2005, Darren.Weber_at_radiology.ucsf.edu
%                    added options to plot surface normals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.4 $ $Date: 2005/08/14 21:22:42 $';
fprintf('FREESURFER_PLOT_SURF [v %s]\n\n',ver(11:15));

if ~exist('normals','var'), normals = 0; end
if isempty(normals), normals = 0; end
if ~exist('orientOnly','var'), orientOnly = 0; end
if isempty(orientOnly), orientOnly = 0; end

if normals,
    facecolor = [0.75 0.1 0.1];
else
    facecolor = [0.65 0.6 0.6];
end

if ~exist('Hf','var'),
    Hf = figure; hold on
end
if isempty(Hf),
    Hf = figure; hold on
end

Hp = patch('vertices',surf.vertices,'faces',surf.faces,...
    'facecolor',facecolor,'edgecolor','none',...
    'FaceLighting','phong');

if normals,
    x = surf.vertices(:,1);
    y = surf.vertices(:,2);
    z = surf.vertices(:,3);
    normals = get(Hp,'VertexNormals');
    if orientOnly,
        [nrm,normals] = colnorm(normals');
        normals = normals';
    end
    u = normals(:,1);
    v = normals(:,2);
    w = normals(:,3);
    quiver3(x,y,z,u,v,w,0) % do not scale the normals
end

daspect([1,1,1])
camlight
camproj('perspective')

material dull
lighting phong

axis off
rotate3d
hold off

return
