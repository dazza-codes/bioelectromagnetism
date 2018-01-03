function emse_write_wfr(file,vertex,face,meshtype,space)

% emse_write_wfr - write mesh to EMSE wireframe (.wfr)
% 
% emse_write_wfr(file,vertex,face,meshtype,space)
% 
% Write a .wfr file, in minor revision 3 format (ascii).
% See the EMSE website at http://www.sourcesignal.com
% for more information on file formats.
% 
% This function assumes the vertex coordinate axes are 
% +X anterior, +Y left, +Z superior
%
% vertex - Nx3 matrix of XYZ values (in meters)
% face - Nx3 matrix of vertex indices for triangulation
% meshtype - a string, with values of:
%
%     'unknown',     
%     'scalp',       
%     'outer skull', 
%     'inner skull', 
%     {'cortex', 'pial', 'white', 'smoothwm'}
%
% space - 'hspace' for head space (electrodes, default)
%         'vspace' for MRI volume space
%


% $Revision: 1.2 $ $Date: 2007/10/29 18:08:47 $


% Copyright (C) 2004  Darren L. Weber
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
%

% History:  12/2004 Darren.Weber_at_radiology.ucsf.edu
%                 - created function from mesh_emse2matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2007/10/29 18:08:47 $';
fprintf('\nEMSE_WRITE_WFR [v%s]\n',ver(11:15));

if ~exist('meshtype', 'var'),
  meshtype = '';
end

if ~exist('space', 'var'),
  space = 'hspace';
end
if isempty(space),
  space = 'hspace';
end


[path,name,ext] = fileparts(file);
ext = '.wfr';
file = fullfile(path,[name ext]);
fprintf('...writing to: %s\n',file);

fid = fopen(file,'w','ieee-le');

if(fid == -1),
  msg = sprintf('...could not open file: %s',file);
  error(msg);
else
  
  % Write prolog
  fprintf(fid,'3\t4000\n');
  fprintf(fid,'3\n');
  
  % Write mesh type
  type = lower(meshtype);
  switch type,
   case 'unknown',
    meshcode = 0;
   case 'scalp',
    meshcode = 40;
   case 'outer skull',
    meshcode = 80;
   case 'inner skull',
    meshcode = 100;
   case {'cortex', 'pial', 'white', 'smoothwm'},
    meshcode = 200;
   otherwise,
    meshcode = 0;
    fprintf('\n...WARNING, unknown meshtype!\n\n');
  end
  
  if strmatch('vspace', space, 'exact'),
    meshcode = meshcode + 80000;
  end
  fprintf(fid, '%d\n', meshcode);
  
  
  % EMSE Voxel Coordinates
  % Voxel coordinates measure location in terms of the voxels inherent in 
  % the given volumetric set. The origin is the bottom (inferior) axial 
  % slice, the posterior row and in the rightmost column. This coordinate 
  % system is right-handed (although, internally, the origin is in the 
  % anterior row, and thus is left-handed; this representation is not 
  % available to the user). The order of the displayed coordinates is 
  % (slice#, row#, column#).
  %
  % EMSE MRI Coordinates
  % MRI coordinates share the same origin as internal voxel coordinates, 
  % but differ from the latter in two ways: first, the coordinates 
  % are measured in millimeters, not voxels. Secondly, the origin is that 
  % of the internal representation; that is, the inferior slice, anterior 
  % row and rightmost column. As mentioned above, this internal representation 
  % is left-handed. To correct for this, the row axis is numbered in the 
  % opposite direction, making the displayed coordinate system right-handed. 
  % The order of the displayed coordinates is (x, y, z).
  
  % Given a point P(x,y,z) in head frame (the activation point on the 
  % cortical mesh) and you want to find the corresponding voxel in the 
  % vmi file.  Symbolically you have P(head) and you want to find P(voxel).
  % 
  % 1.  The registration file contains the matrix HeadToImage,
  %     so P(MRI-mm) = HeadToImage*P(head), where P(MRI-mm) is the 
  %     point in MRI coordinates.
  % 2.  From the voxel size, you can find P(MRI-voxel), which 
  %     is the MRI coordinates expressed in voxels
  % 3.  Use the offset between the MRI coordinate frame and 
  %     the Image coordinate frame to find P(voxel).
  %
  %Demetrios Voreades, Ph.D.
  %Applications Engineer, Source Signal Imaging
  %
  
  
  % Rotate -90 degrees around Z, given that emse coordinates
  % have +X through Nasion and +Y through left ear.
  fprintf('...assuming coordinate axes are +X anterior, +Y left, +Z superior\n');
  %vertex = rz(vertex,-90,'degrees');
  
  % Write vertex data
  for v = 1:size(vertex,1),
    fprintf(fid,'v\t%12.8f\t%12.8f\t%12.8f\n',vertex(v,1),vertex(v,2),vertex(v,3));
  end
  
  % matlab vertex indices start at one,
  % not zero, so we subtract one from matlab values
  fprintf('...subtracting 1 from face indices, so they start at zero\n');
  face = face - 1;
  for t = 1:size(face,1),
    fprintf(fid,'t\t%d\t%d\t%d\t\n',face(t,1),face(t,2),face(t,3));
  end
  
  fclose(fid);
  
end

return
