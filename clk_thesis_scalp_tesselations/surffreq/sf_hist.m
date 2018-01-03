function fh = sf_hist(p, varargin)
% SF_HIST	Plot histograms of surface values
%		SF_HIST(points [, bins] [,fh] [, label])
%		one histogram/column if <= 3 columns
%   	    	  BINS = # of bins in histogram(s)
%   	    	  FH = figure handle to reuse (new figure by default)
%		  LABEL = string for figure title (appended)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - approximate histograms can be produced from sorted input as follows, 
%   though scaling is off by a constant factor
%	sp = sort(p);
%	bsp = sp(linspace(1, size(p, 1), bins + 1));
%	dbsp = diff(bsp);
%	plot(bsp, [(size(p, 1)/bins^2) ./ dbsp ; 0], 'm');


%%% THINGS TO DO
% ? determine # bins from size of input data
% ? deal with multiple array arguments
%   change (S, la) to (S1, la1) and maintain count

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 1) help sf_hist; return; end

%%% interpret optional arguments
bins	= 20;
fh  	= [];
la	= '';

for j = 1:nargin-1
    aj = varargin{j};
    if ischar(aj)		la = aj;
    elseif isempty(aj)
    elseif all(size(aj) == [1 1]) & ishandle(aj) & strcmp(get(aj,'Type'),'figure')
    	fh = aj;
     elseif prod(size(aj)) == 1	bins = aj;
    else
	error(['Unable to process argument.']);
    end
    clear aj;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot histogram(s)

% convert UINT8 to DOUBLE if necessary
if isa(p,'uint8') p = double(p); end

if size(p, 2) > 3, p = p(:); end

if isempty(fh)
    fh = figure( 'Name',	    	'surface hist', ...
    	    	 'PaperPositionMode',	'auto', ...
	    	 'Position',     	[0, 300, size(p, 2)*350, 300]);
    sf_tool axes misc zoom
end

for j = 1:size(p, 2)
    set(0,'CurrentFigure',fh);
    subplot(1, size(p, 2), j);
    hist(p(:, j), bins);
    if (j == 2) title(la); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

