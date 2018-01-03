function entropy = sf_entropy(data, nbins)
% SF_ENTROPY	Compute entropy of data values
%   	    	entropy = SF_ENTROPY(data, nbins)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% H = -   sum    pk   log2(pk  )
% H = - integral p(x) log2(p(x)) dx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

display 'WARNING: SF_ENTROPY has not been tested'

if  (nargin < 2) help sf_entropy; return; end

[hdata,hbins] = hist(data(:), nbins);
hdata = hdata(hdata ~= 0);

entropy = - sum( hdata .* log2(hdata)) * ( mean(diff(hbins)));
