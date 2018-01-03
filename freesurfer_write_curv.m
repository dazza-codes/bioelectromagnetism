function [curv] = freesurfer_write_curv(fname, curv, fnum, ver)

% freesurfer_write_curv - FreeSurfer I/O function to write a curvature file
% 
% [curv] = freesurfer_write_curv(fname, curv, fnum, [ver])
% 
% writes a curvature vector into a binary file
% fname - name of file to write
% curv  - vector of curvatures, one for each vertex
% fnum  - # of faces in surface.
% ver   - 'old' or 'new' version format (new by default)
% 
% See also freesurfer_read_curv, freesurfer_write_surf, freesurfer_write_wfile


if(nargin < 3),
    fprintf('USAGE: curv = freesurfer_write_curv(fname, curv, fnum, [ver])\n');
    return;
end

if ~exist('ver','var'), ver = 'new'; end
if isempty(ver), ver = 'new'; end


% curv array is one value per vertex, 
% so we get the number of vertices:
vnum = length(curv) ;

% open file as a big-endian file
fid = fopen(fname, 'wb', 'b') ;

switch lower(ver),
    
case 'new',
    
    fprintf('...writing new curv version (float)\n');
    
    NEW_VERSION_MAGIC_NUMBER = 16777215;
    freesurfer_fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
    
    fwrite(fid, vnum,'int32');
    fwrite(fid, fnum,'int32');
    fwrite(fid,    1,'int32'); % vals_per_vertex
    fwrite(fid, curv,'float');
    fclose(fid) ;
    
case 'old',
    
    fprintf('...writing old curv version (int16)\n');
    
    freesurfer_fwrite3(fid,vnum) ;
    freesurfer_fwrite3(fid,fnum) ;
    fwrite(fid, curv .* 100, 'int16') ;
end

return
