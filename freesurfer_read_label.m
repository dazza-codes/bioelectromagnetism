function [label] = freesurfer_read_label(sdir, sname, lname)

% label = freesurfer_read_label(sdir, sname, lname)
%
% Reads the label file 'lname' from the subject 'sname' in the subject's
% label directory (given by fullfile(sdir,sname,'label','').  The label
% files for the subject 'sname' are located in the subject's label
% directory.  If 'sdir' is empty, the environment variable 'SUBJECTS_DIR'
% is used.  If they are not defined, sname must contain a full path to the
% subject's files.  'lname' is a cell array of label filenames (with or
% without the .label extension).  If 'sdir' and 'sname' cannot be defined,
% the function attempts to read labels from the current working directory.
%
% Returns the 'label' cell array of structs, as
% label{i}.file - the label file name
% label{i}.data - a matrix, nvertices-by-5, where each column comprises:
% (1)   vertex index into a lh or rh surface, 
% (2-4) xyz at each vertex, given in RAS coordinates for the orig surface,
% (5)   vertex statistic (eg, curvature, thickness, etc.)
%
% For reference, the .label file format is an ascii file:
% The first line is a comment (may have subject name).
% The second line has the number of points (vertices) in the label.
% Subsequent lines have 5 columns:
% 1   - Vertex number (0-based!, this function adds 1 to this vector)
% 2:4 - RAS of the vertex (may be given for the 'orig' surface)
% 5   - Statistic (which can be ignored)
%
% See also: freesurfer_find_labels, freesurfer_label2annotation,
%           freesurfer_read_surf
%

% $Revision: 1.5 $ $Date: 2005/08/25 23:14:06 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2004, Darren.Weber_at_radiology.ucsf.edu
%                    obtained matlab code from the MGH
%           12/2004, Darren.Weber_at_radiology.ucsf.edu
%                    added 1 to vertex index
%           06/2005, Darren.Weber_at_radiology.ucsf.edu
%                    modified help information and the definition of fname
%                    when sname is empty, fname = [lname,'.label'].
%           08/2005, Darren.Weber_at_radiology.ucsf.edu
%                    added sdir input
%                    lname input can be a cell array
%                    label output is a cell array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.5 $';
fprintf('FREESURFER_READ_LABEL [v %s]\n',ver(11:15));

if(nargin < 1)
    help freesurfer_read_label
    return;
end

% ------------------------------------------------
% check inputs and paths

if ~exist('sdir','var'), sdir = []; end
if ~exist('sname','var'), sname = []; end
if ~exist('lname','var'), lname = cell(0); end

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
    msg = sprintf('using pwd, as cannot locate labels path:\n%s', lpath);
    warning(msg);
    lpath = pwd;
end


label = cell(length(lname),1);
for i = 1:length(lname),
    fname = fullfile(lpath,lname{i});
    label{i}.file = fname;
    label{i}.data = read_label(fname);
end

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function label = read_label(fname)

% open it as an ascii file
fid = fopen(fname, 'r');
if(fid == -1)
    msg = sprintf('Could not open %s\n',fname);
    error(msg)
end

fprintf('...reading %s\n',fname);

fgets(fid);  % read the first comment line

% read the number of vertices in the label
line = fgets(fid);
nvertices = sscanf(line, '%d');

fprintf('...reading %d label vertices\n',nvertices);
label = fscanf(fid, '%d %f %f %f %f\n');
fclose(fid);

label = reshape(label, 5, nvertices);
label = label';

fprintf('...adding 1 to vertex index for matlab compatibility\n\n');
label(:,1) = label(:,1) + 1;

return
