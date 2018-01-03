function [pts] = sf_find_3d(data, val, varargin)
% SF_FIND_3D	Find specific values in 3D data set
%		[pts] = SF_FIND_3D(data, val [, options] ...)
%		returns list of occurances of VAL in DAT
%   	    	OPTIONS = {
%   	    	    'edge'  - return coordinates for voxel edges
%   	    	    	      (faster, more points)
%   	    	    'point' - return coordinates for voxel centers
%   	    	    	      (slower, fewer points)
%   	    	    'vol'   - return logical uint8 volume
%   	    	    }

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - rounding boundaries down to nearest voxel
%   - cuts 375k matches to 300k for 4M points
% - acculumating data and postallocating array saves time and memory

%%% THINGS TO DO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 2)
    help sf_find_3d; return
elseif prod(size(val)) ~= 1
    error(['VAL must be a scalar.']);
end

%%% process optional arguments
mode = 'point';

for j = 1:nargin-2
    aj = varargin{j};
    switch aj
    	case {'edge', 'edge2', 'point', 'point2', 'vol'},  mode = aj;
    	otherwise,
    	    error(['Unable to process argument.']);
    end
end
clear aj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare for search
%%% (save data in cells to avoid resizing)

switch mode
    case 'edge',    	tl = 0; td = cell(size(data,3),1);
    case 'point',   	tl = 0; td = cell(size(data,3),1);
    case 'vol',     	pts = logical(uint8(zeros(size(data)))); 
end
d1 = data(:,:,1) >= uint8(val);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% search through dataset    	
for j = 1:size(data,3)
    d0 = d1;
    d1 = data(:,:,j) >= uint8(val);
    
    switch mode
    	case 'edge',
	    [x1,x2] = find(xor(d1(:,[1,1:end-1]  ),d1));
	    [y1,y2] = find(xor(d1(  [1,1:end-1],:),d1));
	    [z1,z2] = find(xor(d0,                 d1));
    	    td{j}   = [ x1 - 0.5, x2,       j(ones(size(x1)));
    		    	y1,       y2 - 0.5, j(ones(size(y1)));
    			z1,       z2,       j(ones(size(z1))) - 0.5];
    	    tl      = tl + size(td{j},1);

    	case 'point',
	    [x1,x2] = find(xor(d1(:,[1,1:end-1]  ),d1) | ...
	    	    	   xor(d1(  [1,1:end-1],:),d1) | ...
	    	    	   xor(d0,                 d1));
	    td{j}   = [x1, x2];
    	    tl      = tl + size(td{j},1);

	case 'vol',
	    pts(:,:,j) = xor(d1(:,[1,1:end-1]  ),d1) | ...
	    	    	 xor(d1(  [1,1:end-1],:),d1) | ...
	    	    	 xor(d0,                 d1);
	otherwise,
	    error(['Unknown mode.']);
    end % switch mode    	
end % for j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% clean up after search
switch mode
    case 'edge',
    	pts = zeros(tl,3);
    	tl = 0;
    	for j=1:size(data,3)
    	    pts(tl+1:tl+length(td{j}),:) = td{j};
    	    tl = tl + length(td{j});
    	end
    case 'point',
    	pts = zeros(tl,3);
    	tl = 0;
    	for j=1:size(data,3)
    	    pts(tl+1:tl+length(td{j}),:) = [td{j} j(ones(length(td{j}),1))];
    	    tl = tl + length(td{j});
    	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
