function freesurfer_write_wfileFast(fname, w)

% freesurfer_write_wfile - FreeSurfer I/O function to write an overlay (*.w) file
%
% [w] = freesurfer_write_wfile(fname, w)
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written, 
%          assumed sorted from vertex 1:N
%
% See also freesurfer_read_wfile, freesurfer_write_surf, freesurfer_write_curv

if(nargin ~= 2)
    msg = sprintf('USAGE: [w] = freesurfer_write_wfile(fname, w)\n');
    error(msg);
end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;

vnum = length(w) ;

fwrite(fid, 0, 'int16') ; % latency integer
freesurfer_fwrite3(fid, vnum) ;   % number of vertices

vertexIndexArray=int32(0:(vnum-1));
vertexBytes=typecast(vertexIndexArray,'uint8');
vertexBytes=reshape(vertexBytes,4,vnum);
%keyboard
dataArray=single(w);
dataBytes=typecast(dataArray,'uint8');
dataBytes=reshape(dataBytes,4,vnum);
dataBytes=flipud(dataBytes);
vertexBytes=vertexBytes(1:3,:);
vertexBytes=flipud(vertexBytes);

outBlock=[vertexBytes;dataBytes];
outBlock=outBlock(:);

fwrite(fid,outBlock,'uint8');

% Replace this loop on ML7 with a typecast version
%
% for i=1:vnum,
%   freesurfer_fwrite3(fid, i-1) ;  % FS vertices start at zero
%   fwrite(fid, w(i), 'float') ;
% end

fclose(fid) ;

return
