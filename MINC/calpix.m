function cp = calpix(x,y)
% CALPIX  Obsolete function - please use PIXELINDEX instead

% $Id: calpix.m,v 1.1 2004/11/12 01:30:50 psdlw Exp $
% $Name:  $

disp ('calpix is obsolete.  Please use pixelindex instead');
cp = pixelindex ([128 128], x, y);
