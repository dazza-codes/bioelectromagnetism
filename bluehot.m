function cmap = bluehot(N);
% BLUEHOT make blue hot to red hot dual intensity colorscale
% function cmap = bluehot(N);
% returns N=64, if N not given.  
% Negative values are blue, positive values are red, with intuitive
%  convention that positive values are "warmer" than negative.  I would have
%  preferred that positive values were coming at you "faster," and therefore
%  had a "blue" shift.  But that's another kind of physics, I guess.

% Copyright (c) 1995, The Regents of the University of California.
% This software was produced under a U.S. Government contract
% (W-7405-ENG-36) by Los Alamos National Laboratory, which is operated
% by the University of California for the U.S. Department of Energy,
% and was funded in part by NIH grant R01-MH53213 through the University
% of Southern California,
% and was funded in part by NIH grant R01-EY08610 to Los Alamos
% National Laboratory.
% The U.S. Government is licensed to use, reproduce, and distribute this
% software.  Permission is granted to the public to copy and use this
% software without charge, provided that this Notice and any statement
% of authorship are reproduced on all copies.  Neither the Government
% nor the University makes any warranty, express or implied, or assumes
% any liability or responsibility for the use of this software.
%
% Author: John C. Mosher, Ph.D.
% Los Alamos National Laboratory
% Group ESA-MT, MS J580
% Los Alamos, NM 87545
% email: mosher@LANL.Gov

% April 4, 1995 author
% April 5, 1995 JCM flipped it for Lewine/Hamalainen convention
% June 8,  1995 JCM, added default N and changed to return size N, not 2N
if(exist('N')~=1),
  N = 64;			% default size
end

cmap1 = hot( ceil(N/2));		% in case its odd
cmap2 = hot(floor(N/2));

cmap = [flipud(cmap1(:,[3 2 1])); cmap2];
% cmap = flipud([flipud(cmap);cmap(:,[3 2 1])]); % near is "hot", neg is "cold"

return
