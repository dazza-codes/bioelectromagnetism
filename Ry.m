function [XYZ] = ry(XYZ,b,units)

% ry - Rotate 3D Cartesian coordinates around the Y axis
%
% Useage: [XYZ] = ry(XYZ,beta,units)
%
% XYZ is a [3,N] or [N,3] matrix of 3D Cartesian coordinates
%
% 'beta'  - angle of rotation about the Y axis
% 'units' - angle is either 'degrees' or 'radians'
%           the default is beta in radians
% 
% If input XYZ = eye(3), the XYZ returned is
% the rotation matrix.
% 
% See also rx rz
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
    b = b*pi/180;
end

Ry = [ cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b) ];

if isequal(size(XYZ,1),3),
    XYZ = Ry * XYZ;
else
    XYZ = XYZ';
    if isequal(size(XYZ,1),3),
        XYZ = [Ry * XYZ]';
    else
        error('Ry: Input XYZ must be [N,3] or [3,N] matrix.\n');
    end
end

return
