% COLORFILE  Defines colors for PGONTRACE
%
% COLOR=COLORFILE(POINT,NORMAL) is an m-file you write
% to define the color of a point on a surface in dependency of its 
% coordinates and its normalvector. (You must, of course, give other
% names than 'Colorfile')
%
% POINT and NORMAL are row vectors containing the points coordinates
% and the components of its normalvector.
% COLOR can be either a 1-by-3 RGB-vector or an index into the
% colormap. In the latter case index need not be scaled.
%
% Example(RGB-vector):
% function color=simplyred(pt,nv);
% color=[1 0 0];
% >>pgontrace(pt,pgon,'simplyred')
%
%
% Example(Index into colormap):
% function color=smooth(pt,nv);
% a=[1;1;2];
% color=abs(nv*a);
% >>pgontrace(pt,pgon,'smooth')
%
% See also  PGONTRACE PGONDISP DOO
