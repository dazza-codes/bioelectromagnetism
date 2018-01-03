function area = sf_gen_area(pts, total, method, p1, p2)
% SF_GEN_AREA	Generate area estimates
%		area = SF_GEN_AREA(pts, total, method [, p1 [, p2]])
%		estimates the area occupied by each value in PTS
%		values are scaled so that they sum to TOTAL
%		  METHOD = {'con', 'edge', 'face', 'rad' }
%		    con:  area = constant for all locations
%		    edge: area = (average_edge_length)^2,   P1 = EDGE list
%		    face: area = (area_of_adjacent_faces),  P1/P2 = FACE list
%		    rad:  area = (min_edge_length)^2,	    

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - edges is a N x 2 matrix of point indices for endpoints
% - faces is a M x 3 matrix of point indices for corners
% - all values are scaled to given total area
%
% - Fechtinger etal comment on various schemes
%   - isolated samples are over-weighted in area-based approach
% - the area-weighting scheme should be independent of extrinsic values, 
%   but probably depends on the number and position of intrinsic values
%
% - if edges correspond to voronoi dual,
%   then face method area equals voronoi area, which is probably best
%   
% - rad2 is alternate radius method which checks distance to all points

%%% THINGS TO DO
% ? implement improved estimates
%   - upper bound on area for isolated points (ala Feichtinger, pg 47)
%     crop to nyquist limit...
% ? adjust/correct area at boundaries

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 3) 
    help sf_gen_area; return;
elseif all(size(pts, 2) ~= [2 3])
    error(['PTS must have 2 or 3 columns, not ', num2str(size(pts,2)), '.']);
elseif ~ischar(method)
    error(['METHOD must be a string.']);
elseif prod(size(total)) ~= 1
    error(['TOTAL must be a scalar.']);
end


switch method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CONSTANT - constant value
    case 'con',
    	area = ones(size(pts, 1), 1);
    	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EDGE - use average length of adjacent edges / 2
    case 'edge',
	if (nargin < 4)
	    error(['no EDGE list.']);
	elseif size(p1, 2) ~= 2
	    error(['EDGE list must have 2 columns, not ',num2str(size(p1, 2)), '.']);
	end

    	epts = reshape(pts(p1,:),[size(p1) size(pts,2)]);
    	%%% ELEN = length of all edges
    	elen = sqrt(sum(diff(epts,1,2).^2,3));
    	area = sf_collect(p1(:), [elen;elen], 'mean', 'full', [size(pts,1) 1]).^2;
    	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FACE - use area of adjacent faces / 3

    %%% triangle area = sqrt(s(s-a)(s-b)(s-c)), 
    %%% where s = (a+b+c)/2, and abc are lengths of sides
    case 'face',
	if (nargin < 4)
	    error(['no FACE list.']);
	elseif size(p1, 2) ~= 3
	    if (nargin < 5)
	    	error(['no FACE list.']);
	    elseif size(p2, 2) == 3
	    	p1 = p2;
	    else 
	    	error(['FACE list must have 3 columns, not ',num2str(size(p1, 2)), '.']);
    	    end
	end

    	fpts = reshape(pts(p1(:,[1 2 3 1]),:),[size(p1,1) 4 size(pts,2)]);
    	%%% FLEN = edge lengths of all faces
    	flen = sqrt(sum(diff(fpts,1,2).^2,3));
    	  fs = sum(flen,2)/2;
    	%%% FA = area of all faces
    	  fa = sqrt(prod([fs, fs(:,[1 1 1])-flen],2));
    	area = sf_collect(p1(:), [fa;fa;fa], 'sum', 'full', [size(pts,1) 1]);
    	
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% RADIUS (radius = 1/2 distance to nearest neighbor)
    % these produce the same results, but the latter are 6X,10X slower
    % (although rad2, rad3 don't require edge information)
    	
    case 'rad',
    	% closest point is self, so we want second closest
    	[ind, area] = sf_find_neighbors(pts,pts,2);
    	area = max(area,[],2);
    	area = area .* area;

    case 'rad1',
	if (nargin < 4)
	    error(['no EDGE list.']);
	elseif size(p1, 2) ~= 2
	    error(['EDGE list must have 2 columns, not ',num2str(size(p1, 2)), '.']);
	end

    	% for each edge, update endpoints if shortest edge so far
	area = Inf; area = area(ones(size(pts, 1), 1));
	for j = p1'
	    area(j) = min(area(j), sum((pts(j(1), :) - pts(j(2), :)).^2));
	end

    case 'rad2',
	area = zeros(size(pts, 1), 1);
	for j = 1:size(pts, 1)
	    jval = pts(j, :);
	    sqdist = sum(((pts - jval(ones(size(pts, 1), 1), :)).^2),2);
	    sqdist(j) = Inf;
	    area(j) = min(sqdist);
	end

    case 'rad3',
    	rpts = pts(:,:,ones(size(pts,1),1));
    	area = squeeze(sum((rpts-permute(rpts,[3 2 1])).^2,2));
    	area(logical(speye(length(area)))) = Inf;
    	area = min(area,[],2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise,
	error(['Unknown METHOD ''', method, '''.']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% scale area to specified total
area = area * total / sum(area);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

