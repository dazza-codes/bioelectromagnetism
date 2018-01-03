function [mri] = mri_open(mri)

% mri_open - function to call various mri data tools
%
% Usage: [mri] = mri_open(mri)
%
% mri is a parameter structure (see mri_toolbox_defaults for
% more details). In this function, it should contain at least
% the following string fields:
%       
%       mri.path - the directory location of the file to load
%       mri.file - the name of the file to load
%       mri.type - the file format (Analyze, FreeSurfer)
%       
%       Analyze is documented in AVW_*_READ
%       FreeSurfer: http://surfer.nmr.mgh.harvard.edu/
%       
% The return structure creates or updates p.mri.data, which contains:
%       
%       mri.data.hdr     struct, eg see avw_hdr_read
%       mri.data.img     3D matrix of image values
%       
% To plot the data returned, set p.mri.plot = 1 before loading, or use:
%       
%       avw_view(mri.data)
%       
% See also, avw_img_read, cor_img_read, avw_view
% 

% $Revision: 1.2 $ $Date: 2008/05/04 18:23:43 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu
%           11/2002, Darren.Weber_at_radiology.ucsf.edu
%                    corrected some bugs and mistakes on p.mri.type
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '$Revision: 1.2 $';
fprintf('MRI_OPEN [v %s]\n',version(11:15));

if ~exist('mri','var'),
  mri = mri_toolbox_defaults;
  fprintf('...creating default p structure.\n');
elseif isempty(mri),
  mri = mri_toolbox_defaults;
  fprintf('...creating default p structure.\n');
end

type = lower(mri.type);

switch type,

case 'analyze',
	
    [path,name,ext] = fileparts(strcat(mri.path, filesep, mri.file));
    file = fullfile(path,[name ext]);
    
    fprintf('...loading Analyze MRI from:\n... %s\n\n',file);
    
    % see avw_img_read for details about orientation
    switch mri.orient
    case 'auto',                mriOrient = '';
    case 'axial unflipped',     mriOrient = 0;
    case 'coronal unflipped',   mriOrient = 1;
    case 'sagittal unflipped',  mriOrient = 2;
    case 'axial flipped',       mriOrient = 3;
    case 'coronal flipped',     mriOrient = 4;
    case 'sagittal flipped',    mriOrient = 5;
    otherwise,                  mriOrient = '';
    end
    
    [mri.data, mri.IEEEMachine ] = avw_img_read(file, mriOrient, mri.IEEEMachine);
	
case 'brainstorm',
    
    fprintf('...BrainStorm not supported yet\n\n');
    return
    %fprintf('...loading BrainStorm data from:\n... %s\n',file);
    
case {'cor','freesurfer'},
	
    % Get Freesurfer data
    [path,name,ext] = fileparts(strcat(mri.path,filesep,mri.file));
    file = fullfile(path,[name ext]);
    
    [ mri.data, mri.IEEEMachine ] = cor_img_read(path, mri.IEEEMachine);
    
otherwise,
    fprintf('...MRI format: %s\n', mri.type);
    fprintf('...Sorry, cannot load this data format at present.\n\n');
    return;
end

if mri.plot, avw_view(mri.data); end

return
