function [retval] = rv(d1, d2);

% RV returnes the residual variance between measured and simulated data
%	[rv] = rv(measured, simulated)
%

% (c) Robert Oostenveld, 1999

retval = mean((d1-d2).^2) ./ mean(d1.^2);

