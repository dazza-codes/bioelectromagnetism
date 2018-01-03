function status = GE_convertEpirecon(firstfile,PFile,outstem, file_type)
%
% status = GE_convertEpirecon(firstfile,PFile,outstem,file_type)
%
%GE_convertEpirecon
%
% Version 1.0
%
% Converts a series of slices reconstructed with epirecon
% saved as EXXXSYYYI1.MR to EXXXSYYYI(Nslices*Nreps).MR
% outstem is the stem used for naming the Analyze files
% file_type = 0 5.X data
% file_type = 1 LX data (default)
%
% status = 1, error
% status = 0, all is well
%
% Souheil J. Inati  
% Dartmouth College
% January 2002
% souheil.inati@dartmouth.edu
%

if (nargin < 3)
        error('Not enough input arguments.')
        return
end

if nargin < 4
  file_type = 1;  % default is LX
end

% Read all the header info from the file
[rdb_hdr,acq_tab,ex_hdr,se_hdr,im_hdr,sizeofpoolheader] = GE_readRawHeader(PFile, file_type);

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

% The image Dimensions
Nx =  im_hdr.dim_X; % X Voxels
Ny =  im_hdr.dim_Y; % Y Voxels
Nz =  im_hdr.slquant;   % Z Voxels
volSize = [Nx Ny Nz];

% Create the template SPM volume structure
V.fname = '';
V.dim = [volSize spm_type('int16')]; % short ints
V.mat = M;
V.pinfo = [1 0 0]';

% Get some stuff from the reconned image header 
[su_hdr_MR,ex_hdr_MR,se_hdr_MR,im_hdr_MR,pix_hdr_MR,im_offset] = GE_readHeader(firstfile);
depth = pix_hdr_MR.img_depth;

% Loop over the volumes until there is no more looping to be done.
for i = 1:im_hdr.nex
  [imageVol] = GE_readVolumeMR(firstfile, volSize, depth, im_offset, i);

  % Create analyze file
  vol_str = sprintf('000%d',i);
  vol_str = vol_str(length(vol_str)-3:length(vol_str));
  outName = [outstem sprintf('_i%s.img',vol_str)];
  V.fname = outName;
  V = spm_create_image(V);

  % Write out the SPM Volume
  % Don't use spm_write_vol because that sets a scale factor.
  for j = 1:Nz
    V = spm_write_plane(V,squeeze(imageVol(:,:,j)),j);
  end

  % Update to the screen
  fprintf('Wrote %s\n',outName);

end

% Done
fprintf('\nConversion Finished. \n');

return
