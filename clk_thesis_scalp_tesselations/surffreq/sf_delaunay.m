function [ed, fa] = sf_delaunay(pts, max)
% SF_DELAUNAY	Generate edges and faces using Delaunay triangulation
%   	    	[edge,face] = SF_DELAUNAY(pts [, max])
%   	    	generates FACEs and EDGEs from Delaunay triangulation
%   	    	ignores potential edges with length > MAX
%   	    	(results are generally similar to sf_gen_edge & sf_gen_face)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if nargin < 1
    help sf_delaunay; return
elseif size(pts,2) ~= 2
    error(['Invalid number of columns (', num2str(size(pts,2)), ') in INTR.']);
end

fa = delaunay(pts(:,1), pts(:,2));
ed = unique(sort([fa(:,1), fa(:,2); ...
    	    	  fa(:,2), fa(:,3); ...
    	    	  fa(:,1), fa(:,3)], 2), 'rows');

if nargin > 1
    epts = reshape(pts(ed             ,:), [size(ed,1) 2 size(pts,2)]);
    elen = sqrt(sum(diff(epts,1,2).^2,3));
    ed = ed(elen < max,:);
    fpts = reshape(pts(fa(:,[1 2 3 1]),:), [size(fa,1) 4 size(pts,2)]);
    flen = sqrt(sum(diff(fpts,1,2).^2,3));
    fa = fa(all(flen < max,2),:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
