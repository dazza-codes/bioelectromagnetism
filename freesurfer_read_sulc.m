function [sulc, fnum] = freesurfer_read_sulc(fname)

% freesurfer_read_sulc - FreeSurfer I/O function to read an sulc file
%
% [sulc, fnum] = freesurfer_read_sulc(fname)
% 
% reads a binary sulc file into a vector
%
% After reading an associated surface, with freesurfer_read_surf, try:
% patch('vertices',vert,'faces',face,...
%       'facecolor','interp','edgecolor','none',...
%       'FaceVertexCData',sulc); light
% 
% See also freesurfer_read_curv, freesurfer_read_surf, freesurfer_read_wfile


% This is just a wrapper, the .sulc format is the same as a .curv format

[sulc, fnum] = freesurfer_read_curv(fname)

return
