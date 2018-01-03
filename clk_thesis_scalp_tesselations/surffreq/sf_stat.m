function [varargout] = sf_stat(varargin)
% SF_STAT	Compute specified statistic(s) across specified signal(s)
%		[val1 ...] = sf_stat(stat1, [...] sig1, sig2 [...])
%	    STAT = {'max', 'min', 'mean', 'median', 'std', 'var'}

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES:

%%% THINGS TO DO:
% ? skip concatenation if only input signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin == 0) help sf_stat; return; end

for nstat = 1:nargin
    if ~ischar(varargin{nstat}) nstat = nstat-1; break; end;
end

nsig  = nargin - nstat;

if (nstat < 1) error(['Need at least one statistic.']); end
if (nsig  < 2) error(['Need at least two signals.']); end

sigsize = size(varargin{nstat+1});
sigs = [];

for j = 1:nsig
    if any(sigsize ~= size(varargin{j+nstat}))
	error(['Signal ', int2str(j), ' has incompatible size.']);
    end
    sigs = [ sigs; varargin{j+nstat}(:)' ];
end

for j = 1:nstat
    varargout{j} = reshape(feval(varargin{j}, sigs), sigsize);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

