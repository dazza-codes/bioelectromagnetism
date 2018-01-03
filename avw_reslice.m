function avw = avw_reslice(avw,dim,pixdim)

% avw_reslice - Linear interpolation (using interp3)
%
% Usage: avw = avw_reslice(avw,dim,pixdim)
%
% avw is the Analyze struct returned by avw_read
% 
% dim is the new volume dimensions (default, dim = [256,256,256])
% pixdim is the new voxel resolution (default, pixdim = [1,1,1] %mm)
%

% $Revision: 1.1 $ $Date: 2005/08/15 21:56:28 $

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

% History:  04/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ver = '$Revision: 1.1 $ $Date: 2005/08/15 21:56:28 $';
fprintf('AVW_RESLICE [v %s]\n',ver(11:15));

if ~exist('dim','var'), dim = [256,256,256]; end
if isempty(dim), dim = [256,256,256]; end

if ~exist('pixdim','var'), pixdim = [1.0,1.0,1.0]; end
if isempty(pixdim), pixdim = [1.0,1.0,1.0]; end

range = dim .* pixdim;

old_dim = double(avw.hdr.dime.dim(2:4));
old_pixdim = double(avw.hdr.dime.pixdim(2:4));
old_range = old_dim .* old_pixdim;

fprintf('...reslicing from: dim = [%d, %d, %d], pixdim = [%6.3f, %6.3f, %6.3f]\n',old_dim,old_pixdim);
fprintf('...reslicing to:   dim = [%d, %d, %d], pixdim = [%6.3f, %6.3f, %6.3f]\n',dim,pixdim);

tic;

% NOTE:
% if we start at zero, we get dim(x)+1 values, one value for each voxel 'edge'
% but this is incompatible with the interp3 function, so we start from the
% first pixdim value below and we get dim(x) values.
%x = 0:old_pixdim(1):old_dim(1)*old_pixdim(1);


% contruct the old meshgrid
x = old_pixdim(1):old_pixdim(1):old_range(1);
y = old_pixdim(2):old_pixdim(2):old_range(2);
z = old_pixdim(3):old_pixdim(3):old_range(3);
[X,Y,Z] = meshgrid(x,y,z);

% contruct the new meshgrid
xi = pixdim(1):pixdim(1):range(1);
yi = pixdim(2):pixdim(2):range(2);
zi = pixdim(3):pixdim(3):range(3);


% NOTE:
% we use EXTRAPVAL below because the range of the original volume may be
% less than the new volume.  When this is the case, the interp3 function
% will not offset the old volume into the center of the new volume.  For
% this to happen, we have to shift the xi/yi/zi values by half the
% difference in the volume size.  That is:

if old_range(1) < range(1),
    xi_shift = ( old_range(1) - range(1) ) / 2;
    xi = xi + xi_shift;
end
if old_range(2) < range(2),
    yi_shift = ( old_range(2) - range(2) ) / 2;
    yi = yi + yi_shift;
end
if old_range(3) < range(3),
    zi_shift = ( old_range(3) - range(3) ) / 2;
    zi = zi + zi_shift;
end
[XI,YI,ZI] = meshgrid(xi,yi,zi);

% run the interpolation
EXTRAPVAL = 0.0;
VI = interp3(X,Y,Z,avw.img,XI,YI,ZI,'linear',EXTRAPVAL);

avw.hdr.dime.dim(2:4) = int16(dim);
avw.hdr.dime.pixdim(2:4) = single(pixdim);
avw.img = VI;

t = toc; fprintf('done (%5.2f sec)\n',t);

return
