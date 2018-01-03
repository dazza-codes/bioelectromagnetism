function [mri] = mri_toolbox_defaults(command,mri);

% mri_toolbox_defaults - Create, read, write mri_toolbox defaults
%
% Useage:   [mri] = mri_toolbox_defaults(command,[mri])
%
% command  =  'create'
%             'read'
%             'write' or 'write_other'
% mri      =  structure returned by create|read command
%             structure to write for write command
%
% All read and write commands are to a matlab .mat file.  The
% 'write_other' command will write out a data specific archive
% that can be accessed from the recent files list.  It will be
% located in the folder where the eeg data is opened.
% 
% Examples:
% mri = mri_toolbox_defaults; % create new defaults
% mri = mri_toolbox_defaults('read'); % read saved defaults
% mri_toolbox_defaults('write',mri);  % write current mri as default
% mri_toolbox_defaults('write_other',mri); % write current mri data
%
% The write command will write to the mri_toolbox
% installation folder.  The read command will first 
% try to find the paramater file in the mri_toolbox
% installation and otherwise recreates the defaults.
%
% Notes:    Handles parameter structure for the 
%           mri_toolbox routines.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Licence:  GNU GPL, no express or implied warranties
% Created:  01/2002 Darren.Weber@flinders.edu.au
% Modified: 02/2002 Darren.Weber@flinders.edu.au
%           - adapted read/write output from .txt format 
%             to .mat format so it will contain data structures.
%             Hence, it is no longer editable text.
%           08/2002 Darren.Weber@flinders.edu.au
%                   added MRI defaults
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mriversion = '$Revision: 1.1 $';
fprintf('\nMRI_TOOLBOX_DEFAULTS [v %s]\n',mriversion(11:15));

if version('-release') < 10,
    msg = printf('matlab release < version 10, mri_toolbox GUIs may fail.\n');
    warning(msg);
end

if ~exist('command','var'), command = 'create'; end

% try to locate the installation path
mriPath = fileparts(which('mri_toolbox'));
if isempty(mriPath),
    msg = sprintf('Cannot find mri_toolbox on the matlab path.\nPlease use the addpath command.\n\n');
    error(msg);
else
    mriPath = strcat(mriPath,filesep);
end


% try to locate the mri_toolbox installation path
mriPath = fileparts(which('avw_view'));
if isempty(mriPath),
    msg = sprintf('Cannot find mri_toolbox on the matlab path.\nPlease install and use the addpath command.\n\n');
    warning(msg);
    mriPath = mriPath;
else
    mriPath = strcat(mriPath,filesep);
end


switch command
case 'create',
    
    fprintf('...creating mri_toolbox defaults\n\n');
    
    mri = [];
    
    % parameters for 'mri_open' & 'avw_view', see these for more information
    mri.path          = strcat(mriPath,'mri_example_data',filesep);
    mri.file          = 'T1.img';  % SPM T1 template
    mri.type          = 'Analyze'; % Analyze files
    mri.orient        = 'auto';
    mri.series        = 1;
    mri.data          = [];
    mri.fiducials     = []; % MRI fiducial points
    mri.IEEEMachine   = 'ieee-be'; % T1 is big endian MRI data
    mri.plot          = 0;
    
    mri.colormap      = gray(255);
    mri.hold          = 0;   % default for GUI hold checkboxes
    
case 'read'
    
    % Look for default file 
    [path,name,ext] = fileparts(strcat(mriPath,'mri_toolbox_defaults.mat'));
    file = fullfile(path,[name ext]);
    if exist(file) == 2,
        fprintf('...reading mri_toolbox defaults from:\n%s\n\n',file);
        load(file);
    else
        fprintf('...cannot locate mri_toolbox defaults - recreating defaults.\n\n');
        mri = mri_toolbox_defaults;
    end
    
    % verify that path to default files exists;
    % if not, create defaults again (I hope this will avoid 
    % new installation problems)
    mriPathExist = exist(mri.path);
    if ~mriPathExist,
        fprintf('...mri.path does not exist - reinitializing defaults.\n');
        mri = mri_toolbox_defaults('create');
        mri_toolbox_defaults('write',mri);
    end
    
case 'write'
    
    if ~exist('mri','var'),
        msg = sprintf('mri argument is essential: mri_toolbox_defaults(''write'',mri)');
        error(msg);
    end
    
    [path,name,ext] = fileparts(strcat(mriPath,'mri_toolbox_defaults.mat'));
    file = fullfile(path,[name ext]);
    
    save(file,'mri');
    
    fprintf('...saved mri_toolbox defaults to:\n%s\n\n',file);
    
case 'write_other'
    
    if ~exist('mri','var'),
        msg = sprintf('mri argument is essential: mri_toolbox_defaults(''write_other'',mri)');
        error(msg);
    end
    
    [path,name,ext] = fileparts(strcat(mri.path, filesep, mri.file));
    ext = '.mat';
    file = fullfile(path,[name ext]);
    
    [newfile,newpath] = uiputfile(file,'Save MRI_TOOLBOX Dataset to File:');
    if (newfile == 0) & (newpath == 0),
        fprintf('...aborting save\n\n');
    else
        [path,name,ext] = fileparts(strcat(newpath, filesep, newfile));
        ext = '.mat';
        file = fullfile(path,[name ext]);
        save(file,'mri');
        fprintf('...saved mri_toolbox parameters to:\n%s\n\n',file);
        recentfiles = mri_toolbox_recent(file);
    end
    
otherwise
    error('...invalid command\n\n');
end

return
