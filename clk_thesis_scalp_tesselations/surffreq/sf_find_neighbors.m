function [ind, dist, cnt, mean] = sf_find_neighbors(src, tgt, cnt, maxdist)
% SF_FIND_NEIGHBORS  	Find neighbors of source points in target points
%   	    	    	[ind, dst, cnt, mean] = SF_FIND_NEIGHBORS(src, tgt [,cnt [,maxdist]])
%   	    	    	return INDices, DISTance, CNT, and MEAN position
%   	    	    	  for up to CNT neighbors within MAXDIST of each source point
%   	    	    	(misses will have dst=maxdist, ind=0
%   	    	    	(work done by C MEX-file)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

