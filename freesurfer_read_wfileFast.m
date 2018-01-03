function [w,v] = freesurfer_read_wfileFast(fname)

% freesurfer_read_wfile - FreeSurfer I/O function to read an overlay (*.w) file
% 
% [w,vert] = freesurfer_read_wfile(fname)
% 
% reads a vector from a binary 'w' file
%	fname - name of file to read from
%	w     - vector of values
%   vert  - vertex indices
% 
% After reading an associated surface, with freesurfer_read_surf, try:
% patch('vertices',vert,'faces',face,...
%       'facecolor','interp','edgecolor','none',...
%       'FaceVertexCData',w); light
% 
% See also freesurfer_write_wfile, freesurfer_read_surf, freesurfer_read_curv

if (nargin ~= 1),
    msg = sprintf('USAGE: [w,v] = freesurfer_read_wfile(fname)\n');
    error(msg);
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0),
    str = sprintf('could not open w file %s.', fname) ;
    error(str) ;
end

fread(fid, 1, 'int16');  % Skip latency int

vnum = freesurfer_fread3(fid) ;  % Number of non-zero values
fprintf('\n%d vertices found\n:',vnum);

v = zeros(vnum,1) ;
w = zeros(vnum,1) ;

% On ML 7 the command 'typecast' allows us to read in the entire data
% structure as chars, then 'typecast' stuff into ints and floats
% See also freesurfer_write_wfileFast

dataBlock=fread(fid,inf,'uint8');
dataBlock=reshape(dataBlock,7,(length(dataBlock)/7));

vertexIndices=uint8(dataBlock(1:3,:));
dataVals=uint8(dataBlock(4:7,:));
vertexIndices=flipud(vertexIndices);
vertexIndices=[vertexIndices;zeros(1,size(vertexIndices,2))];

dataVals=flipud(dataVals);

v=typecast(vertexIndices(:),'int32');
w=typecast(dataVals(:),'single');
w=double(w);


% for i=1:vnum,
%     v(i) = freesurfer_fread3(fid) ;
%     w(i) = fread(fid, 1, 'float') ;
% end

fclose(fid) ;

if nargout > 1,
    fprintf('...adding 1 to vertex indices for matlab compatibility.\n');
    v = v + 1;
end

% w = zeros(max(v),1);
% w(v+1) = w0;

return
