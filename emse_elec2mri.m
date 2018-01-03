function [V] = emse_elec2mri(elec,reg)

% EMSE_ELEC2MRI - Convert point in head frame to mri coordinates
% 
% [vert] = emse_elec2mri(elec,reg)
% 
% elec   - the Nx3 (X,Y,Z) points in a head frame to be converted
% reg    - a structure containing coordinate transform matrices,
%          which is read using emse_open_reg
% 
% Given a point P(x,y,z) in head frame (eg, an activation point on a 
% cortical mesh) this function will find the corresponding voxel in a 
% vmi file.  Symbolically we have P(head) and want to find P(voxel).
% 
% 1.  The registration file contains the matrix HeadToImage,
%     so P(MRI-mm) = P(head)*HeadToImage, where P(MRI-mm) is the 
%     point in MRI coordinates.
% 2.  From the voxel size, you can find P(MRI-voxel), which 
%     is the MRI coordinates expressed in voxels
% 3.  Use the offset between the MRI coordinate frame and 
%     the Image coordinate frame to find P(voxel).
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber@flinders.edu.au
%                    EMSE details thanks to:
%                    Demetrios Voreades, Ph.D.
%                    Applications Engineer, Source Signal Imaging
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('EMSE_ELEC2MRI: In development. Be careful!\n');

% elec input is Nx3 matrix that should be represented 
% in homogenous coordinates:
elec = [ elec ones(size(elec,1),1) ];

% 1.  The registration file contains the matrix HeadToImage,
%     so P(MRI) = HeadToImage*P(head), where P(MRI-mm) is the 
%     point in MRI coordinates.

% However, I've translated HeadToImage, so we now right-multiply,
% which is consistent with a text book account of the subject.

vert = elec * reg.elec2mri;

% reg.elec2mri is a 4x4 matrix, eg:
%
%  -0.9525    0.0452    0.3012         0
%  -0.0522   -0.9985   -0.0154         0
%   0.3000   -0.0304    0.9534         0
%  -0.1295    0.1299    0.0756    1.0000
%
% The first 3x3 cells are the rotations,
% the last row is the translations, and
% the last column is the scale, if any

% In homogeneous coordinates, the last column
% is the scale factor, usually 1
V(:,1) = vert(:,1) ./ vert(:,4);
V(:,2) = vert(:,2) ./ vert(:,4);
V(:,3) = vert(:,3) ./ vert(:,4);

return
