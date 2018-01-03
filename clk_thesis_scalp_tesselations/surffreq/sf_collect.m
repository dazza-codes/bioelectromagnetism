function result = sf_collect(idx,val,varargin)
% SF_COLLECT	Collect terms into matrix
%   	    	result = SF_COLLECT(idx, val [,maxidx] [,options])
%   	    	returns RESULT matrix with each VAL at given IDX
%   	    	options for handling duplicate values 
%   	    	    'max','min','num','sum','first','last','mean'
%   	    	other options: 'full', 'sparse', size of result matrix
%   	    	defaults: 'sum', 'full'
%    	    	(generalization of SPARSE)
%   	    	(work done by C MEX-file)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

