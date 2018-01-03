function vUnit = vector_unit(v, o)

% vector_unit - computes the unit vector of a 3D Cartesian vector
% 
% vUnit = vector_unit(v,[origin])
%
% v - an Mx3 matrix with M row vectors of (x,y,z) components.
% 
% origin - the coordinate system origin (default is [0,0,0]).  This argument
%          is usually a 1x3 vector for a constant origin.  If it is an Mx3
%          matrix, there can be a different origin for each row vector of v.
%
% vUnit  - an Mx3 matrix of unit vectors for each row vector of v, ie:
%
% vUnit = v[xyz] / sqrt((x-xo).^2 + (y-yo).^2 + (z-zo).^2)
% 

% $Revision: 1.3 $ $Date: 2007/12/13 19:20:19 $

% Licence:  GNU GPL, no express or implied warranties
% History:    12/2003, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('o','var'),
  o = [0, 0, 0];
end
if isempty(o),
  o = [0, 0, 0];
end

if ~exist('v','var'),
  error('no input 3D Cartesian vector');
end
if isempty(v),
  error('empty input 3D Cartesian vector');
end

distance = sqrt((v(:,1) - o(:,1)).^2 + (v(:,2) - o(:,2)).^2 + (v(:,3) - o(:,3)).^2);

vUnit = zeros(size(v));
vUnit(:,1) = (v(:,1) - o(:,1)) ./ distance;
vUnit(:,2) = (v(:,2) - o(:,2)) ./ distance;
vUnit(:,3) = (v(:,3) - o(:,3)) ./ distance;

return
