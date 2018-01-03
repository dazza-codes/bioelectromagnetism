function [Hp] = freesurfer_plot_scalar(surf,scalar)

% freesurfer_plot_scalar
%
% [Hp] = freesurfer_plot_scalar(surf,scalar)
%
% returns a handle to the patch object in Hp.
%


% $Revision: 1.1 $ $Date: 2005/08/14 21:24:19 $

% Copyright (C) 2005  Darren L. Weber
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

% History:  09/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $ $Date: 2005/08/14 21:24:19 $';
fprintf('FREESURFER_PLOT_SCALAR [v %s]\n',ver(11:15));


figure; hold on
Hp = patch('vertices',surf.vertices,'faces',surf.faces,...
    'facecolor','interp','edgecolor','none',...
    'FaceVertexCData',scalar);

colormap(jet); colorbar

daspect([1,1,1])
material dull
lighting phong
axis off

rotate3d

return
