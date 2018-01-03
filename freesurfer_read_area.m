function [area, fnum] = freesurfer_read_area(fname)

% freesurfer_read_area - FreeSurfer I/O function to read an area file
%
% [area, fnum] = freesurfer_read_area(fname)
% 
% reads a binary area file into a vector
%
% After reading an associated surface, with freesurfer_read_surf, try:
% patch('vertices',vert,'faces',face,...
%       'facecolor','interp','edgecolor','none',...
%       'FaceVertexCData',area); light
% 
% See also freesurfer_read_curv, freesurfer_read_surf, freesurfer_read_wfile


% This is just a wrapper, the .area format is the same as a .curv format

[area, fnum] = freesurfer_read_curv(fname)

return
