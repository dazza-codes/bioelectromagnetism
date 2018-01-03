function [XYZ] = rz(XYZ,g,units)

% rz - Rotate 3D Cartesian coordinates around the Z axis
%
% Useage: [XYZ] = rz(XYZ,gamma,units)
%
% XYZ is a [3,N] or [N,3] matrix of 3D Cartesian coordinates
%
% 'gamma' - angle of rotation about the Z axis
% 'units' - angle is either 'degrees' or 'radians'
%           the default is gamma in radians
% 
% If input XYZ = eye(3), the XYZ returned is
% the rotation matrix.
% 
% See also rx ry
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  04/2002, Darren.Weber_at_radiology.ucsf.edu
%                    Developed after example 3.1 of
%                    Mathews & Fink (1999), Numerical
%                    Methods Using Matlab. Prentice Hall: NY.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('units','var'), units = 'radians'; end

% convert degrees to radians
if isequal(units,'degrees'),
    g = g*pi/180;
end

Rz = [ cos(g) -sin(g) 0;  sin(g) cos(g) 0;  0 0 1 ];

if isequal(size(XYZ,1),3),
    XYZ = Rz * XYZ;
else
    XYZ = XYZ';
    if isequal(size(XYZ,1),3),
        XYZ = [Rz * XYZ]';
    else
        error('Rz: Input XYZ must be [N,3] or [3,N] matrix.\n');
    end
end

return
