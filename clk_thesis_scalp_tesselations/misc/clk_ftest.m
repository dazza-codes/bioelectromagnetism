function p = clk_ftest(m1,v1,n1,m2,v2,n2)
% CLK_FTEST	F-TEST for significantly different variances
%		p = CLK_FTEST(mean1, var1, num1, mean2, var2, num2)
%   	    	(mean values are ignored)
%
%		NOTE: be sure to square stddev values

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - see Numerical Recipes 13.4 and 6.3

%%% THINGS TO DO
% - test against known cases?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

if (nargin ~= 6)
    help clk_ftest; return;
end

if (v1 > v2)	f = v1/v2; d1 = n1-1; d2 = n2-1;
else		f = v2/v1; d1 = n2-1; d2 = n1-1;
end

p =	betainc(d1/(d1+d2/f), d1/2, d2/2) - ...
	betainc(d2/(d2+d1*f), d2/2, d1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks



