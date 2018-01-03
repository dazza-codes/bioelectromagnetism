function [annots] = freesurfer_read_annotation(fname)

% freesurfer_read_annotation - reads a binary annotation file into a vector
%
% [annots] = freesurfer_read_annotation(fname)
% 
% See also freesurfer_read_surf, freesurfer_read_curv,
%          freesurfer_read_wfile
%

% $Revision: 1.2 $ $Date: 2005/05/18 23:07:42 $

% Licence:  GNU GPL, no implied or express warranties
% History:  07/2004, Darren.Weber_at_radiology.ucsf.edu
%                    adapted from MGH code from Bruce Fischl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $';
fprintf('FREESURFER_READ_ANNOTATION [v %s]\n',ver(11:15));

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b');
if (fid < 0),
	 str = sprintf('could not open annotation file %s.', fname);
	 error(str);
end

fprintf('...reading annotation file: %s\n', fname);
tic;

vnum = fread(fid, 1, 'int32');
tmp = fread(fid, vnum*2, 'int');
annots = tmp(2:2:vnum*2);

fclose(fid);
t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return
