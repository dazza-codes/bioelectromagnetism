function [cintr, cextr] = sf_comp(intr, extr, varargin)
% SF_COMP	Generate composite surface spectrum estimate
%		[fh | cintr, cextr] = sf_comp(intr, extr [,options ...])
%   	    	CINTR, CEXTR = composite intr/extr coordinates (2 outputs)
%   	    	options:
%   	    	  ERANGE = extrinsic axis range
%   	    	  FH     = figure handle (1 output)
%		  EMODE  = { x y z xp yp zp xyz xyzp
%   	    	    	     max mean median min pow std }
%   	    	  IMODE  = { 1d1 1d2 1dm 1dr 1d 2d }
%   	    	  DMODE  = { all bin cum dist }
%   	    	  CMODE  = { cmax cmin cnum csum cfirst clast cmean }
%   	    	emode: extrinsic values - one subplot for each option
%   	    	imode: intrinsic values - one color for each 1d option
%   	    	dmode: displayed values - one line/point set for each option
%   	    	cmode: collection mode for bin'd and cum'd values

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - theoretically, area estimates should be used,
%   but for regular sampling of integer frequencies, Area === 1

%%% THINGS TO DO
% ? use predetermined bins (for DMODE = cum,bin)
% ? use area estimates (for generality; see above)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 2)
    help sf_comp; return;
elseif size(intr, 2) ~= 2
    error(['INTR must have 2 columns, not ', num2str(size(intr,2)), '.']);
elseif size(extr, 2) ~= 3
    error(['EXTR must have 3 columns, not ', num2str(size(extr,2)), '.']);
elseif size(intr, 1) ~= size(extr, 1)
    error(['INTR and EXTR have incompatible sizes.']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process optional arguments
erange = []; fh = []; 
elist = {}; ilist = {}; dlist = {}; cmode = [];
for j=1:nargin-2
    aj = varargin{j};
    if (prod(size(aj)) == 1 & ishandle(aj) & strcmp(get(aj,'Type'),'figure'))
    	fh = aj;
    elseif ischar(aj)
    	aj = lower(varargin{j});
	switch aj
    	    case 'xyz',
    		elist = {elist{:}, 'x','y','z'};
    	    case 'xyzp',
    		elist = {elist{:}, 'xp','yp','zp'};
    	    case {'x','y','z','xp','yp','zp', 'max','mean','median','min','pow','std'},
    		elist = {elist{:}, aj};

    	    case '1d',
    		ilist = {ilist{:}, '1d1','1d2','1dm','1dr'};
    	    case {'1d1','1d2','1dm','1dr','2d'},
    		ilist = {ilist{:}, aj};

    	    case 'all',
    		dlist = {dlist{:}, 'bin','cum','dist'};
    	    case {'bin','cum','dist'},
    		dlist = {dlist{:}, aj};

    	    case {'cmax','cmin','cnum','csum','cfirst','clast','cmean'},
    		cmode = aj(2:end);
    	    otherwise,
    		error(['Invalid optional string argument.']);
	end % switch
    elseif isempty(aj)
    elseif all(size(aj) == [1 2]) erange = aj;
    else error(['Invalid optional argument.']);
    end % if-else
end % for

%%% use defaults if no options specified
if isempty(elist), elist = {'xp','yp','zp'};   	end
if isempty(ilist), ilist = {'1d1','1d2','1dr'}; end
if isempty(dlist), dlist = {'dist','bin','cum'};    	end
if isempty(cmode), cmode = 'sum';   	    	end

if ~isempty(erange) & erange(1) == erange(2)
    erange = erange + 1000 * [-eps eps];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% create figure

if nargout < 2
    if isempty(fh)
	fh = figure( 'Name',  	    	    	'sf_comp', ...
    	    	     'PaperPositionMode', 	'auto', ...
    	    	     'Position',	    	[50 50 300 * length(elist) 250]);
    else
	set(0,'CurrentFigure',fh); clf;
    end
    sf_tool misc zoom
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% process all permutations of options

dflag = 1;
extr = abs(extr);
emin = -1000 * eps;
emax =  1000 * eps;

for i = 1:length(elist), emode = elist{i};

    switch (emode)
	case 'x',	    	    	extr1 = extr(:,1);
	case 'y',	    	    	extr1 = extr(:,2);
	case 'z',	    	    	extr1 = extr(:,3);
	case 'xp',	    	    	extr1 = extr(:,1).^2;
	case 'yp',	    	    	extr1 = extr(:,2).^2;
	case 'zp',	    	    	extr1 = extr(:,3).^2;
	case {'max','min','std'}, 	extr1 = feval(emode,extr,[],2);
	case {'mean','median'},     	extr1 = feval(emode,extr,2);
	case  'pow',	    	    	extr1 = sum(extr.^2,2);
	otherwise,  	error(['Unknown emode: ', emode, '.']);
    end % switch emode
    	
    for j = 1:length(ilist), imode = ilist{j};
    	switch imode
    	    case '1d1', cstr1 = 'b'; intr1 = abs(intr(:,1));
    	    case '1d2', cstr1 = 'g'; intr1 = abs(intr(:,2));
    	    case '1dm', cstr1 = 'm'; intr1 = max(abs(intr),[],2);
    	    case '1dr', cstr1 = 'r'; intr1 = sqrt(sum(intr.^2,2));
    	    case '2d',  dflag = 2;   intr1 = intr;
    	    otherwise, error(['Unknown imode: ', imode, '.']);
    	end % switch imode
    	
    	if dflag == 1
    	    for k=1:length(dlist), dmode = dlist{k};   
    		switch dmode
    	    	    case 'dist', 
    	    	    	[cintr,idx] = sort(intr1);
    	    	    	 cextr = extr1(idx);
    	    	    	 cstr = [cstr1, '.'];
    	    	    case 'bin',  
    	    	    	[cintr,cextr] = sf_comp_bin(intr1, extr1, cmode);
    	    	    	 cstr = [cstr1, '-'];
    	    	    case 'cum',  
    	    	    	[cintr,cextr] = sf_comp_bin(intr1, extr1, cmode);
    	    	    	 cextr = cumsum(cextr); 
    	    	    	 cextr = (cextr(end) - cextr);
    	    	    	 cstr = [cstr1, '-+'];
    		end % switch dmode
    	        
    	        if nargout < 2
    	            %%% plot data for 1D plot
    	    	    set(0,'CurrentFigure',fh);
    	    	    subplot(1,length(elist),i);
    	            hold on; plot(cintr,cextr,cstr); hold off;
    		    if isempty(erange)
    			emin = min(emin, min(cextr));
    			emax = max(emax, max(cextr));
    		    end
    	        end % if nargout < 2
    	    end % for k
    	
    	elseif dflag == 2
    	    if nargout < 2
    	    	%%% create patch for 2D plots
    		[edge,face] = sf_delaunay(intr1);
    	    	set(0,'CurrentFigure',fh);
    	    	subplot(1,length(elist),i);
    		patch('faces',      	    face, ...
    		      'vertices',   	    intr1, ...
    		      'facevertexcdata',    extr1, ...
    		      'edgecolor',  	    'none', ...
    		      'facecolor',  	    'interp');
    		hold on; plot(intr(:,1),intr(:,2), 'k.'); hold off
    		if isempty(erange)
    		    emin = min(emin, min(min(extr1)));
    		    emax = max(emax, max(max(extr1)));
    		end
    	    end % if nargout < 2	
    	    
    	end % if-else
    end % for j
end % for i

%%% clean up figure 
if nargout < 2
    if ~isempty(erange), emin = erange(1); emax = erange(2); end;
    if dflag == 1
    	imin = min(min(intr1)); 
    	imax = max(max(intr1));
    	for i=1:length(elist)
    	    set(0,'CurrentFigure',fh);
    	    subplot(1,length(elist),i);
    	    axis([imin imax emin emax]);
    	    xlabel([elist{i}, ' frequency']); ylabel('amplitude');
    	end
    elseif dflag == 2
    	for i=1:length(elist)
    	    set(0,'CurrentFigure',fh);
    	    subplot(1,length(elist),i);
    	    axis('equal','square');
    	    caxis([emin emax]);
    	    xlabel(['u1']); ylabel(['u2']);
    	    colorbar;
    	    shading interp;
    	end
    	
    	sf_tool axes cmap shade;
    end
end
	
if nargout == 0,
    cintr = []; cextr = [];
elseif nargout == 1,
    cintr = fh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions

function [cin,cex] = sf_comp_bin(in,ex,mode)
    fmin = min(in);
    fmax = max(in);
    fcnt = sqrt(size(in,1))/2;
    cin = linspace(fmin,fmax,fcnt+1)';
    cex = sf_collect(1+floor(fcnt*(in - fmin)/(fmax-fmin)), ...
    	    	     ex, 'full', mode, [fcnt+1 1]);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

