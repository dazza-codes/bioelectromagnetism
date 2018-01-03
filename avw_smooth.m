function avw = avw_smooth(avw,fwhm)

% avw_smooth - Guassian smoothing
%
% Usage: avw = avw_smooth(avw,fwhm)
%
% avw is the Analyze struct returned by avw_read
% fwhm is an odd integer >= 1 or a tuple of odd integers [x y z], which
% determine the voxels to convolve with the Gaussian kernel 
% (default = 3, a 3x3x3 convolution)
%
% see also smooth3
%

% $Revision: 1.2 $ $Date: 2005/08/14 21:18:55 $

% Copyright (C) 2002  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  08/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('fwhm','var'), fwhm = 3; end
if isempty(fwhm), fwhm = 3; end

fprintf('...gaussian smoothing...'); tic;

avw.img = smooth3(avw.img,'gaussian',fwhm);

t = toc; fprintf('done (%5.2f sec)\n',t);

return
