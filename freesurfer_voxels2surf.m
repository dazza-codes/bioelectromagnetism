function [SurfVertices,T] = freesurfer_surf2voxels(VoxVertices)

% freesurfer_voxels2surf - convert voxel into surface coordinates
%
% [SurfVertices,T] = freesurfer_surf2voxels(VoxVertices)
%
% The input 'VoxVertices' are the MRI voxel indices for a Freesurfer MRI
% volume; they are Nx3 coordinates from the sagittal, coronal and axial
% slice planes, in that order. This function converts the MRI voxel indices
% into the FreeSurfer surface RAS coordinates. 
%
% FreeSurfer MRI volumes are 256^3 voxels, 1mm^3 each.
%
% The MRI volume index has an origin at the left, posterior, inferior
% voxel, such that:
% Sagital increases from left to right (+X Right)
% Coronal increases from posterior to anterior (+Y Anterior)
% Axial   increases from inferior to superior (+Z Superior)
%
% The MRI RAS values have an origin at the middle of the volume, in
% approximately voxel 128, 128, 128.  So, any given voxel or RAS coordinate
% can be transformed into the corresponding RAS or voxel index,
% respectively, with the following transformation (T):
%
% T.fsRAS2Vox = [ [-1 0 0 128]' [0 0 -1  128]' [ 0  1 0 128]' [ 0 0 0 1]' ];
% T.fsVox2RAS = [ [-1 0 0 128]' [0 0  1 -128]' [ 0 -1 0 128]' [ 0 0 0 0]' ];
%
% These transform matrices are designed to be right multiplied with a
% matrix of vertex locations, such as:
%
% SurfVertices = VoxVertices * T.fsVox2RAS;
%
% note that VoxVertices is padded out from an [N,3] matrix of vertices (in
% rows) and their XYZ values in columns, into an [N+1,4] matrix, where the 
% bottom row is zeros and the right column is ones (except the last row).
% So SurfVertices(1:N,1:3) gives back just the transformed vertex
% coordinates (this function does this for you!).
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2004, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('FREESURFER_VOXELS2SURF [v %s]\n',ver(11:15)); tic;

%--------------------------------------------------------------------
% Pad out the VoxVertices to Nx4 matrix

Nvertices = size(VoxVertices,1);

right_column = [ ones( Nvertices, 1 ); 0 ]';

VoxVertices = [ [VoxVertices; 0 0 0]  right_column ];

%--------------------------------------------------------------------
% Convert FreeSurfer RAS values into voxel indices

T.fsRAS2Vox = [ [-1 0 0 128]' [0 0 -1  128]' [ 0  1 0 128]' [ 0 0 0 1]' ];
T.fsVox2RAS = [ [-1 0 0 128]' [0 0  1 -128]' [ 0 -1 0 128]' [ 0 0 0 0]' ];

% T.fsVox2RAS =
% 
%     -1     0     0     0
%      0     0    -1     0
%      0     1     0     0
%    128  -128   128     0

% T.fsRAS2Vox =
% 
%     -1     0     0     0
%      0     0     1     0
%      0    -1     0     0
%    128   128   128     1

SurfVertices = VoxVertices * T.fsVox2RAS;
SurfVertices = SurfVertices(1:Nvertices,1:3);

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return










%--------------------------------------------------------------------------
% TESTING CODE

% %------------------------------------------------------------------------
% % Test the conversion from FreeSurfer RAS to Voxel indices

% % FreeSurfer volume_index is: Sag, Axi, Cor (256^3, 1mm^3 voxels)
% fs.volume_index.nas = [128 130 221];
% fs.volume_index.lpa = [203 152 133];
% fs.volume_index.rpa = [ 49 152 133];
% 
% fs.volume_index.mat = [fs.volume_index.nas; fs.volume_index.lpa; fs.volume_index.rpa ];
% 
% fsVox = [ [fs.volume_index.mat; 0 0 0] [1 1 1 0]']

% eg,
% fsVox =
% 
%    128   130   221     1
%    203   152   133     1
%     49   152   133     1
%      0     0     0     0
%
% % FreeSurfer volume_xyz is +X right, +Y anterior, +Z superior; this
% % is the RAS values from the tkmedit viewer of FreeSurfer, where the volume
% % index values are ordered: Sag, Axi, Cor
% fs.volume_xyz.nas = [  0 93  -2];
% fs.volume_xyz.lpa = [-75  5 -24];
% fs.volume_xyz.rpa = [ 79  5 -24];
% 
% fs.volume_xyz.mat = [fs.volume_xyz.nas; fs.volume_xyz.lpa; fs.volume_xyz.rpa ];
% 
% fsRAS = [ [fs.volume_xyz.mat; 0 0 0] [1 1 1 0]']
% 
% T.fsVox2RAS = [ [-1 0 0 128]' [0 0  1 -128]' [ 0 -1 0 128]' [ 0 0 0 0]' ];
% T.fsRAS2Vox = [ [-1 0 0 128]' [0 0 -1  128]' [ 0  1 0 128]' [ 0 0 0 1]' ];
% 
% % T.fsVox2RAS =
% % 
% %     -1     0     0     0
% %      0     0    -1     0
% %      0     1     0     0
% %    128  -128   128     0
% 
% % T.fsRAS2Vox =
% % 
% %     -1     0     0     0
% %      0     0     1     0
% %      0    -1     0     0
% %    128   128   128     1
% 
% fsVox2ras = fsVox * T.fsVox2RAS;
% fsVox2ras = fsVox2ras(1:3,1:3)
% 
% fsRAS2vox = fsRAS * T.fsRAS2Vox;
% fsRAS2vox = fsRAS2vox(1:3,1:3)
