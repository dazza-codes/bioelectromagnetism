function M = sf_plot(pt, varargin)
% SF_PLOT	Plot intrinsic or extrinsic surfaces, generate movie
%		m = SF_PLOT(points [,axisrange] [,{areas|colors}] ...  
%    	    	    	    [,{edges|faces}] [,ah] [,fh] [,label] [,mode] [,mvals])
%		plots POINTS in 2D or 3D as appropriate
%   	    	using AXISRANGE for all axes
%		using circle of specified AREAS
%		or connected by EDGES or FACES (with COLOR) if available
%   	    	  AH/FH = axis/figure handle to reuse (new by default)
%		  LABEL = string for figure title (appended)
%		  MODE = {'abs','mag','pha','magpha','real','imag','realimag'}
%   	    	    'normals' to include face normals
%		  MVALS = [ # frames, degree range, reps ]
%   	    	return MOVIE generated if nargout == 1, or figure handle

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO
% - make axis argument 1x4 or 1x6, add 1x2 caxis argument
% - complains when x-axis range is zero
%   - probably need to adjust axis/view/zoom calls
% ? deal with multiple sets of arguments
% ? specify colors for faces (0.3 and 0.7)
% ? specify color and linestyle for edges ([1 0 0] '-' and [0 1 0] '--')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 1) 
    help sf_plot; return;
elseif all(size(pt, 2) ~= [2 3])
    error(['Invalid number of columns (', num2str(size(pt, 2)), ') in PT.']);
end

%%% interpret optional arguments
ar = [];
co = [];
ed = [];
fa = [];
ah = [];
fh = [];
la = '';
mv = [];
p1 = real(pt); p2 = [];
nrmflg = 0;

for j = 1:nargin-1
    aj = varargin{j};
    if ischar(aj)
    	switch aj
	    case 'abs',		p1 =   abs(pt);	p2 = [];
	    case 'mag',		p1 =   abs(pt);	p2 = [];
	    case 'pha',		p1 = angle(pt);	p2 = [];
	    case 'magpha',	p1 =   abs(pt);	p2 = angle(pt);
	    case 'real',	p1 =  real(pt);	p2 = [];
	    case 'imag',	p1 =  imag(pt);	p2 = [];
	    case 'realimag',	p1 =  real(pt);	p2 = imag(pt);
	    case 'normals', 	nrmflg = 1;
	    otherwise,		la = [la, aj];
	end
    elseif isempty(aj)
    elseif all(size(aj) == [1 1]) & ishandle(aj) & strcmp(get(aj,'Type'),'axes')
    	ah = aj;
    elseif all(size(aj) == [1 1]) & ishandle(aj) & strcmp(get(aj,'Type'),'figure')
    	fh = aj;
    elseif all(size(aj) == [1 1])		co = aj(ones(size(pt, 1), 1));
    elseif all(size(aj) == [1 2])   	    	    	    	ar = aj;
    elseif all(size(aj) == [1 3]) 		        	mv = aj;
    elseif all(size(aj) == [size(pt,1) 1])	    	    	co = aj;
    elseif ~isempty(ed) & all(size(aj) == [size(pt,1) 2])	co = aj;
    elseif ~isempty(fa) & all(size(aj) == [size(pt,1) 3])	co = aj;
    elseif size(aj, 2) == 2			    	    	ed = aj;
    elseif size(aj, 2) == 3			    	    	fa = aj;
    % if arg j is square and same size as points, take diagonal as color/area
    elseif size(aj, 1) == size(aj, 2) & ...
	   size(aj, 1) == size(pt, 1)		co = full(diag(aj));
    else
	error(['Unable to process argument.']);
    end
end
clear aj;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for all plots 

extr = (size(pt, 2) == 3);

if ~isempty(ah)
    %%% reusing axis
    axes(ah);
    
elseif ~isempty(fh)
    %%% reusing figure, so delete axes
    set(0,'CurrentFigure',fh);
    kids = get(fh,'Children');
    kids = kids([strcmp(get(kids, 'Type'), 'axes')]);
    delete(kids);

else
    %%% create new figure
    fh = figure('Name', 			'sf_plot', ...
	    	'PaperPositionMode',		'auto', ...
	    	'Position',			[0, 300, 300, 300], ...
	    	'DefaultLineLineWidth',		0.25, ...
	    	'DefaultLineMarkerSize',	2);

    if extr				view(-45,30); end
					sf_tool axes misc zoom
    if all(size(fa)) | all(size(co))	sf_tool cmap shade; end
    if extr				sf_tool spin; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot faces

if all(size(fa))

    % use indexed or true color as appropriate
    if size(co, 1) == size(pt, 1)
    	switch size(co, 2)
    	    case 2, 	C1 = [co zeros(size(co, 1), 1)];
    	    case {1,3}, C1 = co;
    	    otherwise,	error(['Invalid color values']);
    	end
    	C2 = C1;
    else
    	C1 =  0.3 * ones(size(pt, 1), 1);
    	C2 =  0.7 * ones(size(pt, 1), 1);
    end
    
    % plot appropriate set(s) of vertices and faces
    ax = gca;
    if all(size(p1))
    	patch('faces',      	    	fa, ...
    	      'vertices',   	    	p1, ...
    	      'facevertexcdata',    	C1, ...
    	      'backfacelighting',   	'reverselit', ...
      	      'facecolor',  	    	get(ax,'defaultsurfacefacecolor'), ...
              'edgecolor',  	    	get(ax,'defaultsurfaceedgecolor'));
        hold on; % so both sets will appear
    end
    if all(size(p2))
    	patch('faces',      	    	fa, ...
    	      'vertices',   	    	p2, ...
    	      'facevertexcdata',    	C2, ...
    	      'backfacelighting',   	'reverselit', ...
      	      'facecolor',  	    	get(ax,'defaultsurfacefacecolor'), ...
              'edgecolor',  	    	get(ax,'defaultsurfaceedgecolor'));
    end
    hold off;
    
    % add colorbar for indexed color
    if isempty(ah) & all(size(co) == [size(pt,1) 1])  colorbar; end
    
    if extr & nrmflg
    	fp = reshape(p1(fa,:), [size(fa,1) 3 size(p1,2)]);
    	fc = squeeze(mean(fp,2));
    	fx = -cross(p1(fa(:,1),:) - p1(fa(:,2),:), ...
    	    	    p1(fa(:,1),:) - p1(fa(:,3),:));
    	hold on; quiver3(fc(:,1),fc(:,2),fc(:,3), fx(:,1),fx(:,2),fx(:,3), 'r-');
    end % if nrmflg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot edges

elseif size(ed)
    ax = gca;
    if all(size(p1))
    	patch('faces',	    	    	ed, ...
    	      'vertices',   	    	p1, ...
    	      'edgecolor',  	    	[1 0 0], ...
    	      'linestyle',  	    	'-', ...
    	      'facelighting',	    	'none', ...
    	      'edgelighting',	    	'flat');
    	hold on; % so both sets will appear
    end
    if all(size(p2))
    	patch('faces',	    	    	ed, ...
    	      'vertices',   	    	p2, ...
    	      'edgecolor',  	    	[0 1 0], ...
    	      'linestyle',  	    	'--', ...
    	      'facelighting',	    	'none', ...
    	      'edgelighting',	    	'flat');
    end
    hold off;
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot area circles centered at points

elseif all(size(co))
    co = sqrt(co/pi);
    cpts = [cos(linspace(0, 2*pi, 20)) NaN];
    spts = [sin(linspace(0, 2*pi, 20)) NaN];  
    if all(size(p1))
	xpts = (p1(:, 1*ones(21,1)) + co*spts)';
	ypts = (p1(:, 2*ones(21,1)) + co*cpts)';
	plot(xpts(:), ypts(:), 'r-');
	text(p1(:,1), p1(:,2), ...
	     num2str((1:size(p1,1))'), ...
             'HorizontalAlignment','center');
    	hold on; % so both sets will appear
    end
    if all(size(p2))
	xpts = (p2(:, 1*ones(21,1)) + co*spts)';
	ypts = (p2(:, 2*ones(21,1)) + co*cpts)';
	plot(xpts(:), ypts(:), 'g--');
	text(p2(:,1), p2(:,2), ...
	     num2str((1:size(p2,1))'), ...
             'HorizontalAlignment','center');
    end
    hold off;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot points

else
    if all(size(p1))
    	if extr plot3(p1(:,1), p1(:,2), p1(:,3), 'r*');
    	else	 plot(p1(:,1), p1(:,2),          'r*');
    	end
    	hold on; % so both sets will appear
    end
    
    if all(size(p2))
    	if extr plot3(p2(:,1), p2(:,2), p2(:,3), 'g+');
    	else	 plot(p2(:,1), p2(:,2),          'g+');
    	end
    end
    hold off;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for all plots

axis('equal', 'vis3d');
title(la);

if extr
    if isempty(ar) ar = [ min(axis) max(axis) ]; end
    view(-45,30);
    axis([ ar ar ar ]);
    zoom fill;
    rotate3d on;
else
    if isempty(ar) ar = [ min(axis) max(axis) ]; end
    axis([ ar ar ]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate movie if requested

if all(size(mv))
    num = mv(1); ran = mv(2); rep = mv(3);
    
    axis('off');
    % get size of frame in pixels 
    set(gca, 'Units', 'pixels');
    rect = get(gca, 'Position');
    M = moviein(num, gcf, rect);
    
    for j = 1:num
	view(-45 + j*ran/num, 30);
	M(:,j) = getframe(gcf, rect);
    end
    
    if (nargout == 0) 
    	disp('hit a key to start movie...')
    	pause
    	movie(gcf, M, -rep, 12, rect)
    	clear M
    end
else
    M = fh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

