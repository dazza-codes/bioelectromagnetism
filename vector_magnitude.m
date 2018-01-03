function vMagnitude = vector_magnitude(v,origin)

% vector_magnitude - computes the magnitude of a 3D Cartesian vector
% 
% vMagnitude = vector_magnitude(v,[origin])
%
% v - an Mx3 vector with M rows of (Vx,Vy,Vz) components
% 
% origin - the coordinate system origin, a 1x3 vector.
%          This argument is optional, the default is (0,0,0)
%
% vMagnitude - Mx1 column vector of magnitudes for each row 
%              vector in v, ie:
%
%              sqrt( (Vx-xo).^2 + (Vy-yo).^2 + (Vz-zo).^2 );
%
% The matlab function 'norm' does not give this (as at v6.5).
% 

% $Revision: 1.2 $ $Date: 2007/04/13 18:19:45 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    the function should handle an offset origin
%                    and explicit use of a matrix v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('origin','var'),
  origin = [0, 0, 0];
end
if isempty(origin),
  origin = [0, 0, 0];
end

if ~exist('v','var'),
  error('no input 3D Cartesian vector');
end
if isempty(v),
  error('empty input 3D Cartesian vector');
end

Vx = v(:,1) - origin(1);
Vy = v(:,2) - origin(2);
Vz = v(:,3) - origin(3);

vMagnitude = sqrt( Vx.^2 + Vy.^2 + Vz.^2 );

return
