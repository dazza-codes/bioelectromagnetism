function freesurfer_write_label(sname, lname, label)
% freesurfer_write_label(<sname>, lname, label)
%
% writes the label file 'lname' from the subject 'sname' in the subject's
% label directory from the 'label' matrix, where 'label' will be
% nvertices-by-5 and each column means:
%
% (1) vertex number, (2-4) xyz at each vertex, (5) stat
%
% The Label file format (ascii):
% 
% The first line is a comment (may have subject name).
% The second line has the number of points in the label.
% Subsequent lines have 5 columns:
% 
% 1. Vertex number (0-based!)
% 2-4. RAS of the vertex
% 5. Statistic (which can be ignored)
% 
% This function assumes the 'label' matrix was created with
% freesurfer_read_label, so it will subtract 1 from the first column.
%
% see also, freesurfer_label2annotation, freesurfer_read_surf
%

% $Revision: 1.2 $ $Date: 2004/12/10 19:07:36 $

% Licence:  GNU GPL, no implied or express warranties
% History:  12/2004, Darren.Weber_at_radiology.ucsf.edu
%                    adapted from freesurfer_read_label
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $';
fprintf('\nFREESURFER_WRITE_LABEL [v %s]\n',ver(11:15));

if(nargin ~= 3)
  fprintf('freesurfer_write_label(<sname>, lname, label)\n');
  return;
end

if(~isempty(sname))
  sdir = getenv('SUBJECTS_DIR');
  fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname);
else
  fname = lname;
end

% open it as an ascii file
fid = fopen(fname, 'w');
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fprintf('...writing %s\n',fname);

% write the first comment line
count = fprintf(fid, '#!ascii label, created by freesurfer_write_label.m [v %s]\n',ver(11:15));
if count < 1,
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

% write the number of vertices in the label
nvertices = size(label,1);
fprintf(fid,'%d\n',nvertices);

fprintf('...subtracting 1 from vertex index column of label matrix\n');
label(:,1) = label(:,1) - 1;

fprintf('...writing %d label vertices\n',nvertices);
fprintf(fid,'%8d %12.6f %12.6f %12.6f %12.6f\n',label');
fclose(fid);

return
