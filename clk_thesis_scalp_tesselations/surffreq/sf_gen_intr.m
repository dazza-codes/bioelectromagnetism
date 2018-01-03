function [intr, edge, face] = sf_gen_intr(num, org, range, jit)
% SF_GEN_INTR	Generate set of 2D intrinsic locations
%		[intr,edge,face] = SF_GEN_INTR(num, org [, range [, jit]]) 
%		generates NUM*NUM points in (-RANGE/2 ... +RANGE/2)^2 
%		ORG = { 'dist'  - random with constrained minimum distance
%			'rand'  - random positions
%			'rect'  - rectangular array of positions}
%		    append 'q' to suppress messages, 'v' to display status
%		JIT = variance of gaussian jitter for each point in 'rect'
%			 ( multiplied by range/num )
%		defaults:  RANGE=1.0  JIT=0.01
%		edges and faces generated if specified

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - sample range =      max(intr) - min(intr)
%   sample area  = prod(max(intr) - min(intr))
%
% - all methods:
%   - points should run from -range/2 to range/2 - range/num
%     so that sampling is uniform when wrapped or replicated
% - dist method:
%   - min dist is reduced if num^4 attempts fail
% - dirty triangulation method has problems with regular grids,
%   so generate edges/faces directly, add jitter, or use Delaunay

%%% THINGS TO DO
% ? max edge argument to sf_gen_edge should depend on ORG and JITTER...
% - worry about wrap conditions in dist sampling
% - more methods:
%   - weighted methods to concentrate points in specific area
%     (for specific surface types, or special features)
%   - hexagonal spacing with direct edges & faces (avoid triangulation)
%   - delta-spacing - ala Feichtinger
%   - at each pass, add point with least overlap from set of N candidates
%   - formalize constrained spacing - minimum separation between points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_gen_intr; return; end
if (nargin < 4) jit 	= 0.01; end
if (nargin < 3) range	= 1.00; end
if ~ischar(org) error(['ORG argument must be a string']); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% random organization with minimum distance between points

if strcmp(org(1:4), 'dist')
    intr = zeros(num*num, 2);
    cutoff = range/num;
    if ~strcmp(org, 'distq')
	disp(['Original dist cutoff = ', num2str(cutoff) ]);
    end
    newp = range*(rand(1, 2)*(1-1/num) - 0.5);
    intr(1, :) = newp;

    if strcmp(org, 'distv')
	cpts = cos(linspace(0, 2*pi, 20))' * cutoff;
	spts = sin(linspace(0, 2*pi, 20))' * cutoff;
	figure( 'Name', 	'surface plot', ...
		'Position',	[0, 300, 400, 400], ...
		'Resize', 	'off');
	hold on;
	han = plot(newp(1) + spts, newp(2) + cpts, 'r');
	text(newp(1), newp(2), num2str(1), ...
	     'HorizontalAlignment', 'center', ...
	     'FontSize', [6]);
	axis('equal', 'square');
    	% (PAUSE required to get intermediate display)
	pause(0);
    end
    
    for j = 1:num^2-1
	mindist = 0;
	k = 0;
	while mindist < cutoff
	    newp = range*(rand(1, 2)*(1-1/num) - 0.5);
    	    [ind, mindist] = sf_find_neighbors(newp,intr(1:j,:));
	    k = k+1;
	    % prevent infinite loops
	    if (k == num^2) 
	    	cutoff = cutoff * 0.9; k = 0;
		if ~strcmp(org, 'distq')
		    disp(['Point ', num2str(j), ' mindist cutoff = ', num2str(cutoff) ]);
		end
	    end
	end
	intr(j+1, :) = newp;

	if strcmp(org, 'distv')
	    set(han,'Color','w');
	    han = plot(newp(1) + spts, newp(2) + cpts, 'r');
	    text(newp(1), newp(2), num2str(j+1), ...
		 'HorizontalAlignment', 'center', ...
	     'FontSize', [6]);
	    pause(0);
	end
    end
    if strcmp(org, 'distv') 
	hold off; 
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% random organization

elseif strcmp(org(1:4), 'rand')
    intr = range*(rand(num*num, 2)*(1-1/num) - 0.5);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% rectangular organization

elseif strcmp(org(1:4), 'rect')
    [x, y] = meshgrid(linspace(-range/2, range/2 - range/num, num));
    intr = [reshape(x, num*num, 1), ...
	    reshape(y, num*num, 1)];
    %%% add jitter
    intr = intr + randn(size(intr))*(jit*range/num);
        
    if (jit == 0.0)
    	disp('WARNING: jitter == 0, so edges and faces may have errors');
    end

elseif strcmp(org, 'rect1')
    [x, y] = meshgrid(linspace(-range/2, range/2, num+1));
    intr = [reshape(x(1:num, 1:num), num*num, 1), ...
	    reshape(y(1:num, 1:num), num*num, 1)];
    %%% add jitter
    intr = intr + randn(size(intr))*(jit*range/num);
    
    if (jit == 0.0)
    	disp('WARNING: jitter == 0, so edges and faces may have errors');
    end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error(['Unknown organization (', org, ').']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% print sampling values for both domains

if org(end) ~= 'q'
    disp(['sampling:  num=', num2str(num), ...
		    ' ran=', num2str(range), ...
		    ' max=', num2str(0.5*range), ...
		    ' min=', num2str(range/num)]);
    disp([' dual ->   num=', num2str(num), ...
		    ' ran=', num2str(num/range), ...
		    ' max=', num2str(0.5*num/range), ...
		    ' min=', num2str(1/range)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% generate edges and faces if required

if (nargout > 1)
    [edge,face] = sf_delaunay(intr, 5*range/num);
end
    
if (nargout > 1) & 0
    if org(end) ~= 'q' disp(['generating edges...']); end
    edge = sf_gen_edge(intr, 5*range/num); 
end
if (nargout > 2) & 0
    if org(end) ~= 'q' disp(['generating faces...']); end
    face = sf_gen_face(intr, edge); 
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
