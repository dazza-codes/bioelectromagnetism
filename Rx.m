function [XYZ] = rx(XYZ,a,units)

% rx - Rotate 3D Cartesian coordinates around the X axis
%
% Useage: [XYZ] = rx(XYZ,alpha,units)
%
% XYZ is a [3,N] or [N,3] matrix of 3D Cartesian coordinates
%
% 'alpha' - angle of rotation about the X axis
% 'units' - angle is either 'degrees' or 'radians'
%           the default is alpha in radians
% 
% If input XYZ = eye(3), the XYZ returned is
% the rotation matrix.
% 
% See also ry rz
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
    a = a*pi/180;
end

Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a) ];

if isequal(size(XYZ,1),3),
    XYZ = Rx * XYZ;
else
    XYZ = XYZ';
    if isequal(size(XYZ,1),3),
        XYZ = [Rx * XYZ]';
    else
        error('Rx: Input XYZ must be [N,3] or [3,N] matrix.\n');
    end
end

return
