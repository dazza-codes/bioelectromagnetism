function freesurfer_labels2annotation(sname,hemi,label_file_list,annot_file)

% freesurfer_labels2annotation(sname,hemi,label_file_list,annot_file)
%
% sname = subject name
% hemi = 'lh'  or 'rh'
% label_file_list = a cell array of strings naming the labels, without the
%                   .label suffix
% annot_file = the name of the output annotation file
%
% see also freesurfer_read_label, freesurfer_read_annotation
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2004, Daniel Goldenholz at NMR MGH HARVARD
%           08/2004, Darren Weber, adapted to eeg_toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('FREESURFER_LABELS2ANNOTATION [v %s]\n',ver(11:15));

sdir = getenv('SUBJECTS_DIR');
if isempty(sname),
    sname = getenv('SUBJECT');
end

% step 1: establish number of vertices by loading a surface file
surfname = sprintf('%s/%s/surf/%s.inflated',sdir,sname,hemi);
vertices = freesurfer_read_surf(surfname);
nverts = size(vertices,1);
clear vertices;

% step 2: generate annotation vector
vnums = 0:(nverts-1);
inds = 1:2:nverts*2;
annot_vec = zeros(1,nverts*2);
annot_vec(inds) = vnums;
cnums = zeros(1,nverts);

% step 3: color in the annotation
for i=1:length(label_file_list)
    fname = label_file_list{i} ;
    label_mat = freesurfer_read_label(sname,fname);
    verts = label_mat(:,1);
    thiscolor = ceil(2^24 * i / length(label_file_list)); 
    cnums(verts+1) = thiscolor;
end

% step 4: fill in colors
inds = 2:2:nverts*2;
annot_vec(inds) = cnums;

% step 5: save annotation
% open it as a big-endian file
fid = fopen(annot_file, 'wb', 'b');
if (fid < 0)
    str = sprintf('could not open annotation file %s.', fname);
    error(str);
end
fwrite(fid, nverts, 'int32');
fwrite(fid, annot_vec, 'int');
fclose(fid);

return
