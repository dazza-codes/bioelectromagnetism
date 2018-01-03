function mesh_write_emse(p)

% mesh_write_emse - Save mesh to EMSE (.wfr) file
% 
% USEAGE: mesh_write_emse(p)
% 
% Write a .wfr file, in minor revision 3 format (ascii), for each mesh in
% p.mesh.data.  If any cells of p.mesh.data.meshtype are empty, these cells
% will be skipped.
% 
% See the EMSE website at http://www.sourcesignal.com for more information
% on file formats.
% 
% This function wrapper calls emse_write_wfr.m; it uses the default head
% space coordinates (electrodes).  For more control over coordinates, see
% emse_read_reg.m, emse_elec2mri.m, emse_mri2elec.m, and emse_write_wfr.m
%

% $Revision: 1.3 $ $Date: 2007/10/29 19:35:38 $

% Licence:  GNU GPL, no implied or express warranties
% History:  06/2002 Darren.Weber_at_radiology.ucsf.edu
%                 - created function from mesh_emse2matlab
%           - see cvs logs for all further updates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nMESH_WRITE_EMSE...\n');

if ~exist('p','var'),
  error('...no input p struct.\n');
elseif isempty(p),
  error('...input p struct is empty.\n');
elseif isempty(p.mesh.data),
  error('...input p struct has no mesh data.\n');
end

[path,name,ext] = fileparts(strcat(p.mesh.path,filesep,p.mesh.file));
ext = '.wfr';
file = fullfile(path,[name ext]);

fprintf('...writing EMSE meshes to:\n\t%s\n',fullfile(path,name));

tic;

Meshes = p.mesh.data.meshtype;

for i=1:size(Meshes,2),
  
  if Meshes{i},
    vertex = p.mesh.data.vertices{i};
    face = p.mesh.data.faces{i};
    meshtype = Meshes{i};
    
    % Rotate -90 degrees around Z, given that EMSE electrode coordinates have
    % +X through Nasion and +Y through left ear.
    fprintf('...rotating coordinate axes so +X anterior, +Y left, +Z superior\n');
    vertex = rz(vertex,-90,'degrees');
    
    % create a unique file name for each meshtype
    mesh_file_name = strcat(name,'_',meshtype);
    mesh_file = fullfile(path,[mesh_file_name ext]);
    
    % assume the 'space' option is the default head space
    emse_write_wfr(mesh_file,vertex,face,meshtype)
  end
  
end

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return
