function [Hf,Hp] = freesurfer_plot_curv(surf,curv)

% freesurfer_plot_curv - create a patch surface colored by curvature
%
% [Hf,Hp] = freesurfer_plot_curv(surf,curv)
%
% surf is a struct returned by freesurfer_read_surf, it contains faces and
% vertex coordinates (Nx3, Cartesian xyz).
% curv is an array of Nx1 potential values or Nx3 RGB values
% [Hf,Hp] is a handle to the figure and patch object, respectively
%

% $Revision: 1.2 $ $Date: 2005/07/05 23:36:39 $

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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2005/07/05 23:36:39 $';
fprintf('FREESURFER_PLOT_CURV [v %s]\n',ver(11:15));

if size(curv,2) == 3,
    % we have RGB colors that must be scaled between 0 and 1
    a = 0; b = 1;
    curv = (curv+1)*(b-a)/2+a;
end

%curvZ = (curv - mean(curv)) ./ std(curv)

%posCurvIndex = find(curv >= 0);
%negCurvIndex = setdiff(1:length(curv),posCurvIndex);
%curv(posCurvIndex) = 1;
%curv(negCurvIndex) = -1;

Hf = figure; hold on
Hp = patch('vertices',surf.vertices,'faces',surf.faces,...
    'FaceVertexCData',curv,...
    'facecolor','interp',...
    'edgecolor','none');
daspect([1,1,1])

set(gca,'CLim',[-1.2 1.2]);
colormap(flipud(gray(20)));
colorbar

material dull
lighting phong

axis off
rotate3d

return
