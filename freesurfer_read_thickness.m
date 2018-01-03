function [thickness, fnum] = freesurfer_read_thickness(fname)

% freesurfer_read_thickness - FreeSurfer I/O function to read a thickness file
%
% [thickness, fnum] = freesurfer_read_thickness(fname)
% 
% reads a binary thickness file into a vector
%
% After reading an associated surface, with freesurfer_read_surf, try:
% patch('vertices',vert,'faces',face,...
%       'facecolor','interp','edgecolor','none',...
%       'FaceVertexCData',thickness); light
% 
% See also freesurfer_read_curv, freesurfer_read_surf, freesurfer_read_wfile


% This is just a wrapper, the .thickness format is the same as a .curv format

[thickness, fnum] = freesurfer_read_curv(fname);

return
