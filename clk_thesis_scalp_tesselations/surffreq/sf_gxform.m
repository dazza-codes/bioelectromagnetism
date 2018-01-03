function [opnts, tform] = sf_gxform(ipnts, varargin)
% SF_GXFORM	Apply arbitrary geometric xforms to 2D/3D coordinates
%		outpnts = SF_GXFORM(inpnts [, xform, value] ...)
%		produces OUTPNTS from INPNTS by successive application of
%		XFORM with specified VALUE
%		    XFORM = {'nrm', 	    	    (normalize size and position)	    	
%   	    	    	     're', 'rx', 'ry', 'rz'         	(rotate (in radians))
%			     'se', 'sx', 'sy', 'sz', 'sa'   	(scale)
%			     'te', 'tx', 'ty', 'tz'         	(translate)
%   	    	(append 'r' to perform random xform within range)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - random rotate:	-n to +n
% - random scale:	1/n to n (logarithmic)
% - random translate:	-n to +n

%%% THINGS TO DO
% - how to normalize rotations / align major axes?
%   - moment of inertia

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 1) 
    help sf_gxform; return;
elseif (size(ipnts, 2) < 2) | (size(ipnts, 2) > 4) 
    error(['Invalid number of columns (', num2str(size(ipnts, 2)), ') in INPTS.']);
elseif floor(nargin/2) == ceil(nargin/2)
    error(['Must have odd number of arguments']);
end

tform = eye(4);

%%% process transform operations
for j = 1:2:nargin-1
    xt = varargin{j};
    xv = varargin{j+1};
    xf = eye(4);
    
    switch xt
    
	%%%%%%%%%%%%%%%%%%%%
	%%% normalize
    	case 'nrm',  	xv = max(ipnts) - min(ipnts);
    	    	    	xv(xv < mean(xv)/10000) = 1;
    	    	    	xf(1,1) = 1 / xv(1);
	    	    	xf(2,2) = 1 / xv(2);
	    	    	if size(ipnts,2) == 3
	    	    	    xf(3,3) = 1 / xv(3);
	    	    	end
    	  	    	xf(4,:) = [-mean(ipnts) ./ xv, zeros(3 - size(ipnts,1), 1) 1];

	%%%%%%%%%%%%%%%%%%%%
	%%% rotations
	case 're',  	if length(xv) ~= 3
	    	    	    error(['Argument ',num2str(j+1), ' must contain three values.']);
	    	    	end
		    	xf([2 3],[2 3]) = [ cos(xv(1))  sin(xv(1)); -sin(xv(1))  cos(xv(1))];
	    		xf([1 3],[1 3]) = [ cos(xv(2)) -sin(xv(2));  sin(xv(2))  cos(xv(2))];
	    		xf([1 2],[1 2]) = [ cos(xv(3))  sin(xv(3)); -sin(xv(3))  cos(xv(3))];

	case 'rxr',	xv = xv*(rand*2-1);
			xf([2 3],[2 3]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];
	case 'rx',	xf([2 3],[2 3]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];

	case 'ryr',	xv = xv*(rand*2-1);
			xf([1 3],[1 3]) = [ cos(xv) -sin(xv);  sin(xv)  cos(xv)];
	case 'ry',	xf([1 3],[1 3]) = [ cos(xv) -sin(xv);  sin(xv)  cos(xv)];

	case 'rzr',	xv = xv*(rand*2-1);
			xf([1 2],[1 2]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];
	case 'rz',	xf([1 2],[1 2]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];

	case 'rr',      xv = xv*(rand*2-1);
			xf([2 3],[2 3]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];
			xv = xv*(rand*2-1);
			xf([1 3],[1 3]) = [ cos(xv) -sin(xv);  sin(xv)  cos(xv)];
			xv = xv*(rand*2-1);
			xf([1 2],[1 2]) = [ cos(xv)  sin(xv); -sin(xv)  cos(xv)];


	%%%%%%%%%%%%%%%%%%%%
	%%% scalings
	case 'sar',	rv = exp(log(xv)*(rand*2-1));
			xf(1,1) = rv; xf(2, 2) = rv; xf(3, 3) = rv;
	case 'sa',	xf(1,1) = xv; xf(2, 2) = xv; xf(3, 3) = xv;

	case 'se',  	if length(xv) ~= 3
	    	    	    error(['Argument ',num2str(j+1), ' must contain three values.']);
	    	    	end
	    	    	xf(1,1) = xv(1); xf(2,2) = xv(2); xf(3,3) = xv(3);
	    	    	
	case 'sr',	xf(1,1) = exp(log(xv)*(rand*2-1));
			xf(2,2) = exp(log(xv)*(rand*2-1));
			xf(3,3) = exp(log(xv)*(rand*2-1));
    	    	    	
    	case 'sxn', 	xf(1,1) = 1 / (max(ipnts(:,1)) - min(ipnts(:,1)));
	case 'sxr',	xf(1,1) = exp(log(xv)*(rand*2-1));
	case 'sx',	xf(1,1) = xv;

    	case 'syn', 	xf(2,2) = 1 / (max(ipnts(:,2)) - min(ipnts(:,2)));
	case 'syr',	xf(2,2) = exp(log(xv)*(rand*2-1));
	case 'sy',	xf(2,2) = xv;

    	case 'szn', 	xf(3,3) = 1 / (max(ipnts(:,3)) - min(ipnts(:,3)));
	case 'szr',	xf(3,3) = exp(log(xv)*(rand*2-1));
	case 'sz',	xf(3,3) = xv;


	%%%%%%%%%%%%%%%%%%%%
	%%% translations
	case 'te',  	if length(xv) ~= 3
	    	    	    error(['Argument ',num2str(j+1), ' must contain three values.']);
	    	    	end
	    	    	xf(4,:) = [xv 1];
    	
	case 'tr',	xf(4,:) = [ xv*(rand*2-1) xv*(rand*2-1) xv*(rand*2-1) 1];

	case 'txr',	xf(4,1) = xv*(rand*2-1);
	case 'tx',	xf(4,1) = xv;

	case 'tyr',	xf(4,2) = xv*(rand*2-1);
	case 'ty',	xf(4,2) = xv;

	case 'tzr',	xf(4,3) = xv*(rand*2-1);
	case 'tz',	xf(4,3) = xv;


	%%%%%%%%%%%%%%%%%%%%
	otherwise
	    error(['Invalid XFORM ''', xt, '''.']);
    end
    tform = tform * xf;
end

opnts = [ ipnts, ones(size(ipnts, 1), 4 - size(ipnts, 2)) ] * ...
	tform(:, 1:size(ipnts, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
