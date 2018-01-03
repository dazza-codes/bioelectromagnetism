function odata = sf_sample(idata, rate)
% SF_SAMPLE	Return random sampling of data points
%   	    	[odata] = SF_SAMPLE(idata, rate)
%   	    	if 0 < RATE < 1, RATE = fraction of values to return
%   	    	if 1 < RATE < N, RATE = approx number of values to return

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2)
    help sf_sample; return;
elseif prod(size(rate)) ~= 1
    error(['RATE must be a scalar.']);
end

%%% sample input data
if 0 < rate & rate <= 1
    odata = idata(rand(size(idata,1),1) < rate, :);
elseif 1 < rate & rate < size(idata,1)
    odata = idata(rand(size(idata,1),1) < rate/size(idata,1), :);
elseif rate >= size(idata,1)
    odata = idata;
else
    error(['Invalid RATE.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

