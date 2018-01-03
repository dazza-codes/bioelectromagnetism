%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                    %
% Write SPM mat file given the PFile %
%                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = GE_rawCreateSPMmat(PFileName, file_type)
%function M = GE_rawCreateSPMmat(PFileName, file_type)
%
% Put the appropriate translations and rotations to the M matrix
% given the information in the PFile
%
% file_type = 0 5.X data
% file_type = 1 LX data (default)
%
% S. Inati
% Dartmouth College
% Jan. 2002
%

if nargin < 2
  file_type = 1;  % default is LX
end

% Read all the header info from the file
[rdb_hdr,acq_tab,ex_hdr,se_hdr,im_hdr,sizeofpoolheader] = GE_readRawHeader(PFileName, file_type);

% The conversion from pixels to mm
Dims = diag( [im_hdr.dfov/im_hdr.dim_X, ...
	      im_hdr.dfov_rect/im_hdr.dim_Y, ...   % not sure about this one!!
	      im_hdr.slthick + im_hdr.scanspacing ])

% Compute the coordinate system in the image plane
trhc = acq_tab(1).gw_point1; % Top Left Hand Corner of Image
tlhc = acq_tab(1).gw_point2; % Top Right Hand Corner of Image
brhc = acq_tab(1).gw_point3; % Bottom Right Hand Corner of Image

x = trhc - tlhc; x = x./sqrt(x'*x);   % xhat
y = trhc - brhc; y = y./sqrt(y'*y);   % yhat

% The directional normal to the plane
z = acq_tab(end).gw_point1 - acq_tab(1).gw_point1;
z = z./sqrt(z'*z);% zhat

% Build the M matrix for SPM
% M takes a voxel from the image and gives it a coordinate in mm
% On the scanner the voxels start in the top left hand corner of the
% first image.  In SPM they start in the bottom left hand corner,
% so flip y and set the origin to tlhc.
% NB: The voxel (1,1,1) should have position tlhc
Rot = [x, -y, z];
M = eye(4);
M(1:3,1:3) = Rot  * Dims;
M(1:3,4) = tlhc - Rot*Dims*[1;1;1];

return