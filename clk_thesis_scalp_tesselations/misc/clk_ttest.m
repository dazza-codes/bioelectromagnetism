function p = clk_ttest(m1,v1,n1,m2,v2,n2)
% CLK_TTEST	T-TEST for significantly different means
%		for distributions with unequal variances
%		p = CLK_TTEST(mean1, var1, num1, mean2, var2, num2)
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
    help clk_ttest; return;
end

%%% compute Student's t statistic
t = (m1 - m2) / sqrt((v1/n1) + (v2/n2));

%%% compute degrees of freedom (non-integer)
n = ((v1/n1) + (v2/n2))^2 / ( ((v1/n1)^2/(n1-1)) + ((v2/n2)^2/(n2-1)));

%%% compute probability from incomplete beta function
p = 1 - betainc(n/(n+t*t),n/2,1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks



