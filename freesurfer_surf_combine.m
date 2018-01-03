function [FV3] = freesurfer_surf_combine(FV1, FV2)

% freesurfer_surf_combine - combine lh and rh surfaces
%
% [FV] = freesurfer_surf_combine(FV_lh, FV_rh)
% 
% It is assumed that each FV input is a struct:
% FV.vertices = Nx3
% FV.faces = Nx3
%
% After combining surfaces, try:
% patch('vertices',FV.vertices,'faces',FV.faces,..
%       'facecolor',[1 0 0],'edgecolor','none'); light
% 
% See also freesurfer_read_surf, freesurfer_surf_separate
%


% $Revision: 1.3 $ $Date: 2005/08/14 21:23:08 $

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

% History:  05/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.3 $ $Date: 2005/08/14 21:23:08 $';
fprintf('FREESURFER_SURF_COMBINE [v %s]\n\n',ver(11:15));

FV3.vertices = [FV1.vertices; FV2.vertices];
FV3.faces = [FV1.faces; (FV2.faces + length(FV1.vertices))];

return
