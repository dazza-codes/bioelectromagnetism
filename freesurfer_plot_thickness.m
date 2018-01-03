function [Hp] = freesurfer_plot_thickness(surf,thickness)

% freesurfer_plot_thickness
%
% [Hp] = freesurfer_plot_thickness(surf,thickness)
%
% returns a handle to the patch object in Hp.
%


% $Revision: 1.1 $ $Date: 2005/05/18 23:07:42 $

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

ver = '$Revision: 1.1 $ $Date: 2005/05/18 23:07:42 $';
fprintf('FREESURFER_PLOT_THICKNESS [v %s]\n',ver(11:15));


figure; hold on
Hp = patch('vertices',surf.vertices,'faces',surf.faces,...
    'facecolor','interp','edgecolor','interp',...
    'FaceVertexCData',thickness);

set(gca,'CLim',[0 5]);
colormap(jet); colorbar

daspect([1,1,1])
material dull
lighting phong
axis off

rotate3d

return
