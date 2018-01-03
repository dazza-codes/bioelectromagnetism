function [hspace] = emse_mri2elec(vspace, reg)

% EMSE_MRI2ELEC - Convert mri coordinates to points in head frame
% 
% [hspace] = emse_mri2elec(vspace, reg)
% 
% vspace - a struct with a mesh in MRI volume coordinates (mm)
% vspace.vertices - the Nx3 (X,Y,Z) MRI coordinates to be converted
% vspace.faces    - the Nx3 face connectivity of the mesh
%
% reg - a structure containing coordinate transform matrices,
%       which is returned by emse_read_reg.m
% 
% hspace - a struct like vspace in electrode coordinates (meters)
%
% Given a point P(x,y,z) in MRI volume (eg, an fMRI activation overlayed
% onto a high res T1 volume) this function will find the corresponding
% location in the coordinates of the scalp electrodes (head space).
% Symbolically we have P(voxel) and want to find P(head).
% 
% 1.  Use the offset between the MRI coordinate frame and the MRI volume
% coordinate frame to find P(MRI-voxel).
% 2.  Given P(MRI-voxel) and the voxel size, we can find P(MRI-mm), which is
% the MRI coordinates expressed in mm.
% 3.  The registration file contains the matrix ImageToHeadMatrix, so
% P(head) = P(MRI-mm)*reg.mri2elec, where P(MRI-mm) is the point in MRI
% coordinates, in millimeters.  The values in P(head) are in meters.
% 
% This function performs the last calculation, so all the inputs are assumed
% to be correct.  To load an EMSE wireframe (ie, mesh), see emse_read_wfr.m
% and to load a registration file, see emse_read_reg.m
% 
% See also: EMSE_READ_WFR, EMSE_READ_REG, EMSE_ELEC2MRI
% 

% $Revision: 1.2 $ $Date: 2007/10/27 01:25:48 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber@flinders.edu.au
%                    EMSE details thanks to:
%                    Demetrios Voreades, Ph.D.
%                    Applications Engineer, Source Signal Imaging
%           10/2007, modified code from Justin Ales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2007/10/27 01:25:48 $';
fprintf('EMSE_MRI2ELEC [v %s]\n',ver(11:15));

if(nargin < 1)
    help emse_mri2elec;
    return;
end

if size(vspace.vertices,2) ~= 3,
  error('Input vspace is not an Nx3 matrix')
end

% Convert from millimeter to meter units for EMSE hspace
vs = vspace.vertices / 1000;

% vs is an Nx3 matrix that must be represented in
% homogenous coordinates, so we add ones to the last column
nVertices = size(vs,1);
vs = [vs, ones(nVertices,1)];

% Black-magic: We need to switch to EMSE head space coordinate orientation,
% by taking negative y and re-ordering the axes so we get:
% hspace_x = -1 * vspace-y
% hspace_y = vspace-z
% hspace_z = vspace-x
vs(:,2) = -1 * vs(:,2);
vs = vs(:,[2 3 1 4]);

% Apply the mri2elec transform
hspace.faces = vspace.faces;
hspace.vertices = vs * reg.mri2elec;

% Notes:
% reg.mri2elec is a 4x4 matrix, eg:
%
%  -0.9525    0.0452    0.3012         0
%  -0.0522   -0.9985   -0.0154         0
%   0.3000   -0.0304    0.9534         0
%  -0.1295    0.1299    0.0756    1.0000
%
% The first 3x3 cells are the rotations,
% the last row is the translations,
% the last column is projections (usually 0),
% and the value at 4,4 is the homogenous
% coordinate scale unit, usually 1.

% In homogeneous coordinates, the last column
% is the scale factor, usually 1, but in case
% it is ~= 1
hspace.vertices(:,1) = hspace.vertices(:,1) ./ hspace.vertices(:,4);
hspace.vertices(:,2) = hspace.vertices(:,2) ./ hspace.vertices(:,4);
hspace.vertices(:,3) = hspace.vertices(:,3) ./ hspace.vertices(:,4);

hspace.vertices = hspace.vertices(:,1:3);

return
