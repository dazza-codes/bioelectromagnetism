function x = random_range(a,b,m,n)
% x = random_range(a,b,m,n)
%
% This is wrapper for:
%
% x = a + (b-a) * rand(m,n);
%
% so x is an MxN matrix of random numbers with a 
% uniform distribution on the interval [a,b].
%

% reset the random number generator
rand('state', sum(100*clock));

x = a + (b-a) * rand(m,n);

return
