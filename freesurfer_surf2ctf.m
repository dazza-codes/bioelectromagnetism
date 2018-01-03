function [ctf_shape,T] = freesurfer_surf2ctf(SurfVertices)

% freesurfer_surf2ctf - convert a freesurfer surface into CTF shape
%
% [ctf_shape,T] = freesurfer_surf2ctf(vertices)
%
% The input 'vertices' are obtained from freesurfer_read_surf; they are Nx3
% vertex coordinates in the freesurfer RAS coordinate framework. This
% function first converts the vertex coordinates into the FreeSurfer MRI
% voxel index coordinate system and then from there into the CTF MRI voxel
% coordinate system.  This assumes that a freesurfer MRI volume has been
% loaded into CTF MRI format correctly.
%
% To output the results as a CTF .shape file, use:
% ctf_write_headshape(ctf_shape,'filename.shape');
%
% Both FreeSurfer and CTF volumes are 256^3 voxels, 1mm^3 each.  However,
% they are different orientations, so the same voxel in each format has
% different voxel indices.
%
% CTF MRI default volume index has an origin at the left, anterior, superior
% voxel, such that:
% Sag increases from left to right (+X Right)
% Cor increases from anterior to posterior (+Y Posterior)
% Axi increases from superior to inferior (+Z Inferior)
%
% The origin is at the left, anterior, superior voxel
% CTF volume_index is: Sag, Cor, Axi (256^3, 1mm^3 voxels)
% CTF Sag = 256 - FreeSurfer Sag
% CTF Axi = FreeSurfer Axi
% CTF Cor = 256 - FreeSurfer Cor
%
% FreeSurfer volume_index is: Sag, Axi, Cor (256^3, 1mm^3 voxels)
% FreeSurfer Sag = 256 - CTF Sag
% FreeSurfer Axi = CTF Axi
% FreeSurfer Cor = 256 - CTF Cor
%
% So, we can convert from one to the other with the following
% transformations (T):
%
% T.fsVox2ctf = [ -1 0 0 0; 0 0  1 0; 0 -1 0 0; 256 256   0 1 ];
% T.ctf2fsVox = [ -1 0 0 0; 0 0 -1 0; 0  1 0 0; 256   0 256 1 ];
%
% This function calls freesurfer_surf2voxels, which employs the transforms:
%
% T.fsRAS2Vox = [ [-1 0 0 128]' [0 0 -1  128]' [ 0  1 0 128]' [ 0 0 0 1]' ];
% T.fsVox2RAS = [ [-1 0 0 128]' [0 0  1 -128]' [ 0 -1 0 128]' [ 0 0 0 0]' ];
%
% All transforms are designed to be right multiplied, like this:
%
% ctfVox = fsRAS * T.fsRAS2Vox * T.fsVox2ctf;
%
% see also freesurfer_surf2voxels, ctf_mri2head
%

% $Revision: 1.2 $ $Date: 2005/06/27 21:24:36 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%--------------------------------------------------------------------
% Convert FreeSurfer Surface Vertices from RAS to voxel indices
VoxVertices = freesurfer_surf2voxels(SurfVertices);



ver = '$Revision: 1.2 $';
fprintf('FREESURFER_SURF2CTF [v %s]\n',ver(11:15)); tic;

%--------------------------------------------------------------------
% Pad out the VoxVertices to Nx4 matrix

Nvertices = size(VoxVertices,1);

right_column = [ ones( Nvertices, 1 ); 0 ];

VoxVertices = [ [VoxVertices; 0 0 0]  right_column ];

%--------------------------------------------------------------------
% Convert FreeSurfer voxel indices into CTF voxel indices

T.fsVox2ctf = [  -1   0   0   0;  0   0   1   0;  0  -1   0   0; 256  256    0   1 ];
T.ctf2fsVox = [  -1   0   0   0;  0   0  -1   0;  0   1   0   0; 256    0  256   1 ];

ctf_shape = VoxVertices * T.fsVox2ctf;

ctf_shape = ctf_shape(1:Nvertices,1:3);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return












%--------------------------------------------------------------------------
% TESTING CODE


% CTF MRI default volume index has an origin at the left, anterior, superior
% voxel, such that Sag slices increase from left to right, Cor
% slices increase from anterior to posterior, and Axi slices
% increase from superior to inferior

% CTF volume is 256^3 voxels, 1mm^3 each
% Sag increases from left to right (+X Right)
% Cor increases from anterior to posterior (+Y Posterior)
% Axi increases from superior to inferior (+Z Inferior)
% The origin is at the left, anterior, superior voxel
% CTF volume_index is: Sag, Cor, Axi (256^3, 1mm^3 voxels)
% CTF Sag = 256 - FreeSurfer Sag
% CTF Axi = FreeSurfer Axi
% CTF Cor = 256 - FreeSurfer Cor
ctf.volume_index.nas = [128  35 130]; % voxels
ctf.volume_index.lpa = [ 53 123 152];
ctf.volume_index.rpa = [207 123 152];

ctf.volume_index.mat = [ctf.volume_index.nas; ctf.volume_index.lpa; ctf.volume_index.rpa ];

ctfVox = [ [ctf.volume_index.mat; 0 0 0] [ 1 1 1 0 ]' ];

% eg
% ctfmat =
% 
%    128    35   130     1
%     53   123   152     1
%    207   123   152     1
%      0     0     0     0

% FreeSurfer volume_index is: Sag, Axi, Cor (256^3, 1mm^3 voxels)
% FreeSurfer Sag = 256 - CTF Sag
% FreeSurfer Axi = CTF Axi
% FreeSurfer Cor = 256 - CTF Cor
fs.volume_index.nas = [128 130 221];
fs.volume_index.lpa = [203 152 133];
fs.volume_index.rpa = [ 49 152 133];

fs.volume_index.mat = [fs.volume_index.nas; fs.volume_index.lpa; fs.volume_index.rpa ];

fsVox = [ [fs.volume_index.mat; 0 0 0] [1 1 1 0]']

% eg,
% fsVox =
% 
%    128   130   221     1
%    203   152   133     1
%     49   152   133     1
%      0     0     0     0



% This is the transformation matrix from FreeSurfer voxel indices 
% into CTF MRI voxel coordinates.
T.fs2ctf = [  -1     0     0     0  ;
               0     0     1     0  ;
               0    -1     0     0  ;
             256   256     0     1  ];

% This is the transformation matrix from CTF voxel indices into FreeSurfer
% voxel indices
T.ctf2fs = [  -1     0     0     0 ;
               0     0    -1     0 ;
               0     1     0     0 ;
             256     0   256     1 ];

% or just, T.ctf2fs = inv(T.fs2ctf);

fs2ctf = fsVox  * T.fs2ctf;

ctf2fs = ctfVox * T.ctf2fs;



% these .sph values are not the same as ctf.mri.hdr.headOrigin_*
% values; these ones correspond to
% ctf.mri.hdr.HeadModel_Info.defaultSphereXYZ, which are actually
% given in mm, whereas the values below in ctf.volume_xyz.sph are
% in cm (taken from the MRIViewer window
ctf.volume_index.sph = [130 135 103];

% CTF volume_xyz is: +X nasion, +Y left, +Z superior
% convert cm to mm
ctf.volume_xyz.nas = [18.58  0.00 0.00] * 10;
ctf.volume_xyz.lpa = [ 4.38  9.14 0.00] * 10;
ctf.volume_xyz.rpa = [-0.17 -7.70 0.00] * 10;

ctf.volume_xyz.sph = [ 0.00  0.00 5.00] * 10;

% NOTE
% 02/22/04 Actually do not REALLY know what these values relate to!
% The CTF coordinates have different distances between the
% fiducials from the freesurfer coordinates, despite having
% identical voxel dimensions between freesurfer and CTF!
