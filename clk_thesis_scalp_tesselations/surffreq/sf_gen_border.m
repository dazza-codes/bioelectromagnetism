function brdr = sf_gen_border(pnts, method)
% SF_GEN_BORDER	Identify extreme points in dataset
%		brdr = SF_GEN_BORDER(pnts, method)
%		set global SF_GEN_BORDER_BPNTS to border thickness (default == 1.5)
%
%	METHOD =
%	    ALL:	all points are border points
%	    NONE:	no points are border points
%   	    HULL:   	use convex hull to set border
%	    DIS_RAD:	use distance to set radial border
%	    DIS_RECT:	use distance to set rectangular border
%	    NUM_RAD:	use numbers  to set radial border
%	    NUM_RECT:	use numbers  to set rectangular border

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - to plot border points (brdr) or center points (~brdr):
%   sf_plot(pnts([brdr], :))
%   (see end of file)

% - open surfaces: rescaling doesn't work - corners warp more than sides
% - closed surfaces: have coordinate values wrap

%%% THINGS TO DO
% ? select boundary points whose face angles don't total 360


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 2) 
    help sf_gen_border; return;
elseif all(size(pnts, 2) ~= [2 3])
    error(['PNTS must have 2 (or 3) columns,  not ', num2str(size(pnts, 2)), '.']);
end

bpnts = 1.5;

global SF_GEN_BORDER_BPNTS
if prod(size(SF_GEN_BORDER_BPNTS)) == 1  bpnts = SF_GEN_BORDER_BPNTS; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% border methods

switch method
    case 'all',
    	brdr = logical( ones(size(pnts, 1), 1));

    case 'none',
    	brdr = logical(zeros(size(pnts, 1), 1));

    case {'hull','hull1'},
    	brdr = logical(zeros(size(pnts, 1), 1));
    	ch = convhull(pnts(:,1),pnts(:,2));
    	brdr(ch) = 1;

    case {'hull2','hull3','hull4','hull5','hull6','hull7','hull8','hull9','hull10'},
    	brdr = logical(zeros(size(pnts, 1), 1));
    	mpnt = mean(pnts);
    	for j = 1:str2num(method(5:end))
    	    ch = convhull(pnts(:,1), pnts(:,2));
    	    brdr(ch) = 1;
    	    % move CH points to middle so new CH can be found
    	    pnts(ch,:) = mpnt(ones(size(ch)),:) + 20*eps*randn(size(ch,1),2);  	    
    	end
    		    	
    case 'dis_rad',
	%%% radial distance
	bf = bpnts / sqrt(size(pnts, 1));
	imin = min(pnts);
	imax = max(pnts);
	imid = (imin + imax)/2;
	iran = (1-2*bf) * min(imax - imin)/2;
	brdr = sum(((pnts - imid(ones(size(pnts, 1), 1), :)).^2) > iran^2, 2);

    case 'dis_rect',
	%%% rectangular distance
	bf = bpnts / sqrt(size(pnts, 1));
	imin = min(pnts);
	imax = max(pnts);
	iran = (imax - imin) * bf;
	imin = imin + iran;
	imax = imax - iran;
	brdr = pnts < imin(ones(size(pnts, 1), 1), :) | ...
	       pnts > imax(ones(size(pnts, 1), 1), :);
	brdr = any(brdr, 2);
    
    case 'num_rad',
	%%% radial number (round up to integer for indexing)
	bf = ceil((sqrt(size(pnts, 1))/2 - bpnts)^2 * pi);
	imin = min(pnts);
	imax = max(pnts);
	imid = (imin + imax)/2;
	srt = sort(sum(((pnts - imid(ones(size(pnts, 1), 1), :)).^2), 2));
	iran = srt(bf);
	brdr =     sum(((pnts - imid(ones(size(pnts, 1), 1), :)).^2) > iran, 2);

    case 'num_rect',
	%%% rectangular number (round up to integer for indexing)
	bf = ceil(bpnts * sqrt(size(pnts, 1)));
	srt = sort(pnts);
	imin = srt(bf, :);
	imax = srt(size(srt, 1) - bf, :);
	brdr = pnts < imin(ones(size(pnts, 1), 1), :) | ...
	       pnts > imax(ones(size(pnts, 1), 1), :);
	brdr = any(brdr, 2);

    otherwise,
    	error(['Invalid METHOD ''', method, '''.']);
end

% disp([num2str(length(brdr)), ' points, ', num2str(nnz(brdr)), ' on border ']);

if nargout < 1
    figure( 'Name', 	'border plot', ...
	    'Position',	[0, 300, 400, 400]);
    X1 = pnts([ brdr], 1); 
    Y1 = pnts([ brdr], 2);
    X2 = pnts([~brdr], 1);
    Y2 = pnts([~brdr], 2);
    if size(pnts, 2) == 2
	plot(X1, Y1, 'r*', X2, Y2, 'g+');
	axis('equal', 'square');
	sf_tool axes misc zoom
    else
	Z1 = pnts([ brdr], 3);
	Z2 = pnts([~brdr], 3);
	plot3(X1, Y1, Z1, 'r*', X2, Y2, Z2, 'g+');
	axis('equal', 'square');
	zoom fill
	sf_tool axes misc spin
    end
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

