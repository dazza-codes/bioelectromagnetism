function M = sf_splot(intr, extr, varargin)
% SF_SPLOT	Plot surface, generate movie
%		m = SF_SPLOT(intr, extr [,fh] [,interp] [,mode] [,mvals])
%		animates INTR/EXTR surface using
%   	    	  FH = figure handle to reuse (new figure by default)
%   	    	  INTERP = {'linear','cubic','nearest','invdist'}
%   	    	    interpolation method for griddata
%		  CMODE = {'x','y','z','l'}
%   	    	  MODE  = {'abs','mag','pha','real','imag'}
%		  MVALS = [ # frames, degree range, repetitions ]
%   	    	return MOVIE generated, or figure handle

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_splot; return; end
if size(intr, 2) ~= 2
    error(['Invalid number of columns (', num2str(size(intr, 2)), ') in INTR.']);
elseif size(extr, 2) ~= 3
    error(['Invalid number of columns (', num2str(size(extr, 2)), ') in EXTR.']);
elseif size(intr, 1) ~= size(extr, 1)
    error(['INTR and EXTR do not have equal number of rows (', ...
           num2str(size(intr, 1)), ' ~= ', num2str(size(extr, 1)), ').']);
end

%%% interpret optional arguments
cmode = 'l';
fh = [];
interp = 'invdist';
mode = 'real';
mv = [];

for j = 1:nargin-2
    aj = varargin{j};
    if ischar(aj)
    	switch aj
    	    case {'x','y','z','l'}, 
    	    	cmode = aj;
    	    case {'linear', 'cubic', 'nearest', 'invdist'}, 
    	    	interp = aj;
    	    case {'abs','mag','pha','real','imag'},
    	    	mode = aj;
    	    otherwise, 
    	    	error(['Invalid string argument (', aj, '.']);
    	end
    elseif isempty(aj)
    elseif ishandle(aj) & strcmp(get(aj,'Type'),'figure')
    	fh = aj;
    elseif all(size(aj) == [1 3])   
    	mv = aj;
    end
end
clear aj;

%%% compute scaling values
numc = 2.0 *sqrt(length(intr));
[gu,gv] = meshgrid(linspace(min(min(intr)), max(max(intr)),numc));

if isempty(fh)
    fh = figure( 'Name', 		'sf_splot', ...
	    	 'PaperPosition',	[0,   3,   5    5], ...
	    	 'Position', 	    	[0, 300, 500, 500], ...
	    	 'Resize', 	    	'off');
    view(-45,30); sf_tool axes cmap misc shade spin
else
    set(0,'CurrentFigure',fh);
end

switch mode
    case 'abs',     extr = abs(extr);
    case 'mag',     extr = abs(extr);
    case 'pha',     extr = angle(extr);
    case 'real',    extr = real(extr);
    case 'imag',    extr = imag(extr);
    otherwise,      error(['Unknown mode.']);
end

x = griddata(intr(:,1),intr(:,2), extr(:,1), gu,gv, interp);
y = griddata(intr(:,1),intr(:,2), extr(:,2), gu,gv, interp);
z = griddata(intr(:,1),intr(:,2), extr(:,3), gu,gv, interp);

switch cmode
    case 'x', 	    surf(x, y, z, x);
    case 'y',	    surf(x, y, z, y);
    case 'z', 	    surf(x, y, z, z);
    otherwise,	    surfl(x, y, z);
end

view(-45, 30); 
axis('equal', 'square', 'vis3d'); 
zoom fill;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate movie if requested

if size(mv)
    num = mv(1); ran = mv(2); rep = mv(3);

    axis('off');
    % get size of frame in pixels 
    set(gca, 'Units', 'pixels');
    rect = get(gca, 'Position');
    M = moviein(num, gcf, rect);

    for j = 1:num
	view(-45 + ran*j/num, 30);
	M(:,j) = getframe(gcf, rect);
    end

    if (nargout == 0)
    	disp('hit a key to start movie...')
    	pause
    	movie(gcf, M, -rep, 24, rect)
    	clear M
    end
else
    M = fh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

