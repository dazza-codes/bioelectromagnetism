function center = avw_center(avw)

% avw_center - find center of a volume
%
% [center] = avw_center(avw)
%
% avw  - an Analyze 7.5 data struct, see avw_read
%
% center is the Cartesian coordinates for
% the volume center.  It is a struct with
% coordinates in both voxels and mm (1x3).
% 
% center.corner is the result of floor(xdim/2),
% center.abs is the result of (xdim/2).
% 
% Finding the center of a voxel based
% volume can be done in several ways.
% 
% The absolute center point in a volume will 
% lie either within the center of the middle
% voxel or at the boundary between two voxels
% (depending on whether the volume has an odd
% or even number of voxels in any dimension).
% 
% The corner values are the voxel that lies 
% just before the absolute center point of 
% the volume.
% 
% The mm coordinates are simply the voxel values
% multiplied by the pixel dimensions.
%
% Given correctly oriented Analyze 7.5 files, the 
% corner values lie at the right, posterior 
% and inferior corner of the voxel.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.1 $]';
fprintf('\nAVW_CENTER [v%s]\n',version(12:16));  tic;

% Extract info from avw.hdr
xdim = double(avw.hdr.dime.dim(2));
ydim = double(avw.hdr.dime.dim(3));
zdim = double(avw.hdr.dime.dim(4));

xpix = double(avw.hdr.dime.pixdim(2));
ypix = double(avw.hdr.dime.pixdim(3));
zpix = double(avw.hdr.dime.pixdim(4));

% Find center voxel of volume
center.corner.voxels = zeros(1,3);
center.corner.voxels(1,1) = floor(xdim/2);
center.corner.voxels(1,2) = floor(ydim/2);
center.corner.voxels(1,3) = floor(zdim/2);

% Find mm coordinates of that voxel
center.corner.mm = zeros(1,3);
center.corner.mm(1,1) = center.corner.voxels(1,1) .* xpix;
center.corner.mm(1,2) = center.corner.voxels(1,2) .* ypix;
center.corner.mm(1,3) = center.corner.voxels(1,3) .* zpix;

% Find center voxel of volume
center.abs.voxels = zeros(1,3);
center.abs.voxels(1,1) = xdim/2;
center.abs.voxels(1,2) = ydim/2;
center.abs.voxels(1,3) = zdim/2;

% Find mm coordinates of that voxel
center.abs.mm = zeros(1,3);
center.abs.mm(1,1) = center.abs.voxels(1,1) .* xpix;
center.abs.mm(1,2) = center.abs.voxels(1,2) .* ypix;
center.abs.mm(1,3) = center.abs.voxels(1,3) .* zpix;

return
