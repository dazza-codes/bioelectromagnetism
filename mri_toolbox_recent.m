function [ files, varargout ] = mri_toolbox_recent(filename,command);

% mri_toolbox_recent - Keep track of mri_toolbox .mat files
% 
% Useage:   [files,mri] = mri_toolbox_recent([filename],[command])
% 
% filename = path & filename of an mri_toolbox .mat file to
%            be added to the recent files list
% command  = 'load', to load one of the recent files
%            'clear', to clear all recent files
% 
% If filename is not given, the function returns
% the most recent files saved (if any).  If it 
% is given, the filename is added to the list.
% 
% If the load command is given, filename.mat is
% loaded and the mri struct is returned.  Otherwise
% the mri struct is empty.
% 
% Tracks recent mri structure saves to .mat files
% for the mri_toolbox gui.  The return variable
% 'files' is a cell array of strings containing
% the filesystem location of recent mri files saved.
% 
% The recent files list is limited to 20.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Licence:  GNU GPL, no express or implied warranties
% Created:  07/2003 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LIMIT = 20;

if ~exist('filename','var'),
    filename = '';
elseif isempty(filename),
    filename = '';
end
if ~exist('command','var'),
    command = '';
elseif isempty(command),
    command = '';
end


% -- locate the installation path and the recent file list

path = fileparts(which('mri_toolbox_recent'));
if isempty(path),
    msg = sprintf('Cannot find mri_toolbox on the matlab path.\n\nPlease use the addpath command.\n\n');
    error(msg);
else
    path = strcat(path,filesep);
end


% -- load the recent 'files' cell array

recentfile = strcat(path,'mri_toolbox_recent.mat');
if exist(recentfile) == 2,
    load(recentfile);
else
    files = cell(1);
end


% -- add a new file to the 'files' array

if ~isempty(filename),
    if isempty(files{1}),
        files{1} = filename;
    else
        % should check to see if already in list
        % if so, rearrange the list
        
        i = strmatch(filename,char(files),'exact');
        if ~isempty(i),
            % remove it from current list
            if isequal(i,1),         files = files(2:end);
            elseif isequal(i,LIMIT), files = files(1:end-1);
            else                     files = files([1:i-1 i+1:end]);
            end
        end
        files = fliplr(files);   % put most recent at end, for now
        files{end+1} = filename; % add filename to end
        files = fliplr(files);   % now put most recent at beginning
        
    end
    
    % only keep 'LIMIT' most recent
    if size(files,2) > LIMIT,
        files = files(1:LIMIT);
    end
    
    save(recentfile,'files');
end


switch command,
case 'load',
    if ~isempty(filename),
        mri = [];
        load(filename);
    end
case 'clear',
    files = cell(1);
    save(recentfile,'files');
otherwise,
end



if (nargout > 1) & exist('mri','var'),
    if ~isempty(p),
        varargout{1} = mri;
    end
end



return
