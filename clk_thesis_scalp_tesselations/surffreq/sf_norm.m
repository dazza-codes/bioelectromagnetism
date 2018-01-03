function val = sf_norm(points, type, varargin)
% SF_NORM	Compute per-column (1,2,fro) norms if <= 3 columns,
%		regular matrix norms if > 3 columns
%		val = SF_NORM(points, type [, area] [, exp] [, label] )
%		  TYPE 'energy' computes total energy in points
%		  TYPE 'power' computes average power in points
%		  TYPE '1', '2', 'inf', 'fro' passed to MATLAB norm
%		weighted by sample AREA if available
%		raised to EXP power (default == 2) for energy/power
%		results printed using LABEL (appended) if present

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - norm requires full matrices for some norm types

%%% THINGS TO DO
% ? accept and return multiple norm types 
% ? how should area values interact with norms other than energy/power

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

if     (nargin < 2) help sf_norm; return;
end

area 	= ones(size(points, 1), 1);
exp	= 2;
label 	= '';

%%% process optional arguments
for j = 1:nargin-2
    pn = varargin{j};
    if ischar(pn)				label = [label, pn];
    elseif all(size(pn) == [1 1])		  exp = pn;
    elseif all(size(pn) == [size(points, 1) 1])	 area = pn;
    end
end


%%% matrix norms if # cols > 3
if (size(points,2) > 3)
    if     strcmp(type, 'energy')
    	val = sum(sum(area(:, ones(size(points,2),1)) .* abs(points).^exp));
    elseif strcmp(type, 'power')
    	val = sum(sum(area(:, ones(size(points,2),1)) .* abs(points).^exp)) ...
    	      ./ prod(size(points));
    else
    	val = norm(full(points), type);
    end


%%% per-column norms if # cols <= 3
else
    if     strcmp(type, 'energy')
    	val = sum(area(:, ones(size(points,2),1)) .* abs(points).^exp);
    elseif strcmp(type, 'power')
    	val = sum(area(:, ones(size(points,2),1)) .* abs(points).^exp) ...
    	      ./ size(points,1);
    else
	val = zeros(1, size(points, 2));
	for j = 1:size(points, 2)
	    val(j) = norm(full(points(:, j)), type);
	end
    end
end

%%% print label
if ~isempty(label)
    fprintf(1, '%s-norm(%s) = %8.5f %8.5f %8.5f', type, label, val);
    fprintf(1, '\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

