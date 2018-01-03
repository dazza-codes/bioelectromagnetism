function [labels] = freesurfer_find_labels(sdir, sname, lname)

% labels = freesurfer_find_labels([sdir], [sname], [lname])
%
% Finds all the label files for the subject 'sname' in the subject's
% label directory.  If 'sdir' is empty, the environment variable
% 'SUBJECTS_DIR' is used.  If they are not defined, sname may contain a
% full path to the subject's files.  Otherwise, 'sname' is the directory
% name of a subject in the freesurfer 'SUBJECTS_DIR', which contains a
% subdirectory called 'label', which is searched for all files that contain
% the .label file extension.  If all of the above is undefined, the current
% working directory (pwd) is searched for all *.label files.  The optional
% argument 'lname' is a substring of the label filenames, for selection of
% a subset of the label files.
%
% Returns 'labels', a cell array of label file names.  Each of these
% filenames can be input to freesurfer_read_label.
%
% See also: freesurfer_read_label, freesurfer_label2annotation
%

% $Revision: 1.1 $ $Date: 2005/08/25 23:14:06 $

% Copyright (C) 2005  Darren L. Weber
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

% Licence:  GNU GPL, no implied or express warranties
% History:  09/2005, Darren.Weber_at_radiology.ucsf.edu
%                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('FREESURFER_FIND_LABELS [v %s]\n',ver(11:15));

% ------------------------------------------------
% check inputs and paths

if ~exist('sdir','var'), sdir = []; end
if ~exist('sname','var'), sname = []; end
if ~exist('lname','var'), lname = []; end

if isempty(sdir),
    sdir = getenv('SUBJECTS_DIR');
end
if isempty(sdir),
    [sdir,sname,ext] = fileparts(sname);
end
if ~exist(sdir,'file'),
    warning('cannot locate sdir, using current working directory');
    sdir = pwd;
end

lpath = fullfile(sdir,sname,'label','');

if ~exist(lpath,'file'),
    msg = sprintf('using pwd, labels path does not exist:\n%s', lpath);
    warning(msg);
    lpath = pwd;
end


% ------------------------------------------------
% find label files

direc = dir(lpath);

labels = cell(0);
filenames = cell(0);
[filenames{1:length(direc),1}] = deal(direc.name);

% check for the .label extension
Index = zeros(length(filenames),1);
for i = 1:length(filenames),
    if strfind(filenames{i},'.label'),
        Index(i) = 1;
    end
end
filenames = filenames(find(Index));

if lname,
    Index = zeros(length(filenames),1);
    for i = 1:length(filenames),
        % screen for 'lname' files
        if strfind(filenames{i},lname),
            Index(i) = 1;
        end
    end
    filenames = filenames(find(Index));
end

labels = filenames;

if isempty(labels),
    msg = sprintf('cannot locate any labels');
    warning(msg);
end

return
