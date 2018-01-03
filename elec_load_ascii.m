function [elec,X,Y,Z,theta,phi,r] = elec_load_ascii(filename,coordinates,format,xo,yo,zo,N)

% elec_load_ascii - Read ascii electrode file
% 
% Load an ascii electrode file and return the electrode labels and
% coordinates.  This function can read and return Cartesian (x,y,z) and/or
% spherical (theta,phi,r) coordinates.  The only difference between this
% function and elec_load is that this function does not require 'type' in
% the data file and there is some flexibility to define the data format
% specifier in this function.  An electrode label is assumed, to provide a
% facility to define different electrode types.  If you have no electrode
% labels, add a constant text identifier to the start of every line in the
% input data file.
%
% The ascii file format is defined by the 'format' input.  It is assumed
% that each row of the file comprises an electrode label, followed by the
% 'Cartesian' (x,y,z) coordinates (in cm) or the spherical (theta,phi,r) in
% degrees or radians ('Spherical1' and 'Sperical2' coordinates).  Each field
% is separated by white space.  It is assumed that all rows of the input
% file contain only scalp electrode coordinates (there is no differentiation
% of electrode types - such as fiducials, EOG, reference, etc. - apart from
% their label).  For example,
%
% Fz Xcm Ycm Zcm
% 
% Usage:
% 
% [elec,X,Y,Z,theta,phi,r] = elec_load_ascii(filename,coordinates,format,xo,yo,zo,N)
%
% where:
% 
% filename - 'path\filename'
%
% (xo,yo,zo) - the origin, (0,0,0) by default.  If one electrode label is
% 'origin', the values of that "electrode" will be used instead.
%
% N - how many electrodes to read (128 by default).
%
% coordinates - a string option for input data type:
%
%   'Cartesian'  = cartesian in cm (this is default);
%                  format is ['label' X Y Z].
%   'Spherical1' = spherical, theta and phi in degrees,
%   'Spherical2' = spherical, theta and phi in radians,
%                  format is ['label' theta phi r]
%
%  theta - counterclockwise rotation from +x in x-y plane (azimuth),
%  phi   - elevation with respect to z-axis,
%  r     - radius in cm,
%  +ve x-axis from origin through T4 (right ear)
%  +ve y-axis from origin through Nasion (theta = 90 deg = pi/2 rad)
%
% format - a format string to read the data.  Each line of the file should
% contain values separated by white space.  The usual c specifications apply
% (see help textread for more detail).  For example, a file that contains
% rows of cartesian values like this:
%
% Fz 0.0001 8.6728 5.8761
%
% would be read with a format specification as follows:
%
% [elec,X,Y,Z] = textread(file,'%s %f %f %f', N);
%
% If the values were spherical, then we have:
%
% [elec,theta,phi,r] = textread(file, '%s %f %f %f', N);
%
% So, this function is just a fancy wrapper for textread!
%
% Result:
%
% elec is an electrode label, as cellstr array (Nx1)
% (X,Y,Z) & (theta,phi,r) as double floating point (Nx1), with
% theta and phi in radians.  The origin is (0,0,0), unless forced
% otherwise with the input arguments.
%
% Notes:
% 
% ii) Conversion from spherical to Cartesian coordinates is:
% 
% x = r .* sin(phi) .* cos(theta);
% y = r .* sin(phi) .* sin(theta);
% z = r .* cos(phi);
% 
% Phi here is elevation with respect to z-axis.  The matlab
% function sph2cart uses elevation with respect to xy plane.
%
% iii) Conversion from Cartesian to spherical coordinates is:
%
% theta = atan2( y, x );
% phi = atan2( sqrt( x.^2 + y.^2 ), z );
% r = sqrt( x.^2 + y.^2 + z.^2);
% 
% Phi here is the elevation with respect to z-axis.  The matlab
% function cart2sph uses elevation with respect to xy plane.
%

% $Revision: 1.1 $ $Date: 2007/03/15 17:27:16 $

% Licence:  GNU GPL, no express or implied warranties
% History:  03/2007, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('N', 'var'),
  N = 128;
end
if isempty(N),
  N = 128;
end
if ~exist('coordinates', 'var'),
  coordinates = 'cartesian';
end
if isempty(coordinates),
  coordinates = 'cartesian';
end
if ~exist('format','var'),
  format = '%s %f %f %f';
end
if isempty(format),
  format = '%s %f %f %f';
end
if ~exist('xo', 'var'),
  xo = 0;
end
if isempty(xo),
  xo = 0;
end
if ~exist('yo', 'var'),
  yo = 0;
end
if isempty(yo),
  yo = 0;
end
if ~exist('zo', 'var'),
  zo = 0;
end
if isempty(zo),
  zo = 0;
end

tic;

  [path,name,ext] = fileparts(filename);
  file = fullfile(path,[name ext]);

  eegversion = '$Revision: 1.1 $';
  fprintf('ELEC_LOAD_ASCII [v %s]\n',eegversion(11:15));
  fprintf('...loading ''%s'' electrodes from:\n\t%s\n', coordinates, file);

  switch lower(coordinates),
    
   case 'cartesian',
    
    [elec,X,Y,Z] = textread(file, format, N);
    
    fprintf('...converting from cm to meters.\n');
    X = X ./ 100;
    Y = Y ./ 100;
    Z = Z ./ 100;
    
    % If possible, adjust all coordinates so that origin is (0,0,0)
    index = find(elec == 'origin');
    if ~isempty(index),
      xo = X(index);
      yo = Y(index);
      zo = Z(index);
    end
    fprintf('...centering origin at (0,0,0).\n');
    X = X - xo;
    Y = Y - yo;
    Z = Z - zo;
    xo = X(index);
    yo = Y(index);
    zo = Z(index);
    
    % Convert to spherical
    fprintf('...calculating spherical coordinates.\n');
    theta = atan2( (Y-yo), (X-xo) );
    phi = atan2( sqrt( (X-xo).^2 + (Y-yo).^2 ), (Z-zo) );
    r = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2);
    
   case 'spherical1',                   % degrees
    
    [elec,theta,phi,r] = textread(file, format, N);
    % convert theta and phi to radians
    theta = theta * (pi/180);
    phi = phi * (pi/180);
    [X, Y, Z] = elec_sph2cart(theta,phi,r,1);
    fprintf('...converting from cm to meters.\n');
    X = X ./ 100;
    Y = Y ./ 100;
    Z = Z ./ 100;
    
   case 'spherical2',                   % radians
    
    [elec,theta,phi,r] = textread(file, format, N);
    [X, Y, Z] = elec_sph2cart(theta,phi,r,0);
    fprintf('...converting from cm to meters.\n');
    X = X ./ 100;
    Y = Y ./ 100;
    Z = Z ./ 100;
    
   otherwise
    
    doc elec_load;
    msg = sprintf('...invalid coordinate type: ''%s'', see ''help elec_load''\n', coordinates);
    error(msg);
    
  end

  fprintf('...loaded %d electrodes\n', size(X,1));

  t = toc;
  fprintf('...done (%6.2f sec).\n\n',t);

  return
