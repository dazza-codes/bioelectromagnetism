function label = freesurfer_labels2TalDaemon(label)

% label = freesurfer_labels2TalDaemon(label)
%
% 'label' is a cell array of structs returned from
% freesurfer_read_label.  This function parses all the label{i}.data
% through the Talairach Daemon (given an appropriate installation) and
% returns the results into label{i}.td (a struct with fields of data
% returned by the Talairach Daemon).
%
% See also: freesurfer_find_labels, freesurfer_read_label,
%           freesurfer_label2annotation, freesurfer_read_surf
%
% See online: http://ric.uthscsa.edu/TDinfo
%             http://www.brainmap.org/
%

% $Revision: 1.1 $ $Date: 2005/08/25 23:14:06 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2005, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('FREESURFER_LABELS2TALDAEMON [v %s]\n',ver(11:15));

if(nargin < 1)
    help freesurfer_labels2TalDaemon
    return;
end

% identify subject and find the talairach matrix
[labelPath,labelFile,ext] = fileparts(label{1}.file);
subjPath = strrep(labelPath,[filesep,'label'],'');
TALxfmFile = fullfile(subjPath,'mri','transforms','talairach.xfm');
TALxfm = freesurfer_read_talxfm(TALxfmFile);

for i = 1:length(label),
    fname = label{i}.file;
    RASdata = label{i}.data(:,[2:4]);
    [MNIv,TALv] = freesurfer_surf2tal(RASdata,TALxfm);
    fname = strrep(fname,'.label','.tal');
    TalDaemon_write(fname,TALv)
    TDfname = TalDaemon_run(fname);
    label{i}.td = TalDaemon_read(TDfname);
end

return
