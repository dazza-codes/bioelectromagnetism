function extro = sf_filter(intri, extri, varargin)
% SF_FILTER   	Filter extrinsic coordinates
%   	    	extro = SF_FILT(intri, extri [, filter [,p1 ...]])
%   	    	apply sequence of FILTERs (with parameters P1, ...)
%   	    	to given INTRinsic and EXTRinsic Input coordinates
%   	    	and return EXTRinsic Output coordinates
%   	FILTER = {
%   	    box_high_rad    - radial highpass boxcar, BW=P1
%   	    box_high_rect   - rect   highpass boxcar, BW=P1
%   	    box_low_rad     - radial  lowpass boxcar, BW=P1
%   	    box_low_rect    - rect    lowpass boxcar, BW=P1
%   	    but_low_rad     - radial  lowpass Butterworth, BW=p1, order=p2
%   	    but_low_rect    - rect    lowpass Butterworth, BW=p1, order=p2
%   	    }

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
%

%%% THINGS TO DO
% ? edge/face based smoothing (see sf_opt_intr)
% ? other filter windows (see Opp&Sch pg 239ff)
% ? prune filtered points, edges, and faces (use NaN instead of 0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 2)
    help sf_filter; return;
elseif size(intri, 2) ~= 2
    error(['INTRI must have 2 columns, not ', num2str(size(intri, 2))]);
elseif size(extri, 2) ~= 3
    error(['EXTRI must have 3 columns, not ', num2str(size(extri, 2))]);
elseif size(intri, 1) ~= size(extri, 1)
    error(['INTRI and EXTRI must have equal number of rows (', ...
	   num2str(size(intri, 1)), ' ~= ', num2str(size(extri, 1)), ').']);
end


extro = extri;

%%% step through arguments and perform appropriate actions
j=1;
while j <= nargin-2
    aj = varargin{j};
    if ~ischar(aj)
    	error(['Variable argument ', num2str(j), ' is not a string.']);
    end
    switch aj
    	%%% radial highpass boxcar filter
    	case 'box_high_rad',
    	    if (j+1 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    j = j+1;
    	    if ~isnumeric(p1)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    elseif length(p1) == 1
    	    	filt = ((intri(:, 1).^2 + intri(:, 2).^2) > p1^2);
    	    	extro = extro .* filt(:, [1 1 1]);
    	    else
    	    	error(['Invalid parameter for ', aj, '.']);
    	    end

    	%%% rectangular highpass boxcar filter
     	case 'box_high_rect',
    	    if (j+1 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    j = j+1;
    	    if ~isnumeric(p1)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    elseif length(p1) == 1
    	    	filt = (abs(intri(:, 1)) > p1 | ...
    	    	    	abs(intri(:, 2)) > p1);
    	    	extro = extro .* filt(:, [1 1 1]);
    	    elseif length(p1) == 2
    	    	filt = (abs(intri(:, 1)) > p1(1) | ...
    	    	    	abs(intri(:, 2)) > p1(2));
    	    	extro = extro .* filt(:, [1 1 1]);
    	    else
    	    	error(['Invalid parameter for ', aj, '.']);
    	    end

    	%%% radial lowpass boxcar filter
    	case 'box_low_rad',
    	    if (j+1 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    j = j+1;
    	    if ~isnumeric(p1)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    elseif length(p1) == 1
    	    	filt = ((intri(:, 1).^2 + intri(:, 2).^2) < p1^2);
    	    	extro = extro .* filt(:, [1 1 1]);
    	    else
    	    	error(['Invalid parameter for ', aj, '.']);
    	    end

   	%%% rectangular lowpass boxcar filter
   	case 'box_low_rect',
    	    if (j+1 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    j = j+1;
    	    if ~isnumeric(p1)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    elseif length(p1) == 1
    	    	filt = (abs(intri(:, 1)) < p1 & ...
    	    	    	abs(intri(:, 2)) < p1);
    	    	extro = extro .* filt(:, [1 1 1]);
    	    elseif length(p1) == 2
    	    	filt = (abs(intri(:, 1)) < p1(1) & ...
    	    	    	abs(intri(:, 2)) < p1(2));
    	    	extro = extro .* filt(:, [1 1 1]);
    	    else
    	    	error(['Invalid parameter for ', aj, '.']);
    	    end

    	%%% radial lowpass Butterworth filter
    	case 'but_low_rad',
    	    if (j+2 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    p2 = varargin{j+2};
    	    j = j+2;
    	    if ~isnumeric(p1) | ~isnumeric(p2)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    else
	    	filt = sf_butter((intri(:, 1).^2 + intri(:, 2).^2).^(1/2), p1, p2);
	    	extro = extro .* filt(:, [1 1 1]);
	    end
	    
    	%%% rectangular lowpass Butterworth filter
    	case 'but_low_rect',
    	    if (j+2 > nargin)
    	    	error(['Missing parameter(s).']);
    	    end
    	    p1 = varargin{j+1};
    	    p2 = varargin{j+2};
    	    j = j+2;
    	    if ~isnumeric(p1) | ~isnumeric(p2)
    	    	error(['Non-numeric parameter for ', aj, '.']);
    	    elseif length(p1) == 1
	    	filt = sf_butter(intri(:, 1), p1, p2) .* ...
	    	       sf_butter(intri(:, 2), p1, p2);
	    	extro = extro .* filt(:, [1 1 1]);
    	    elseif length(p1) == 2
	    	filt = sf_butter(intri(:, 1), p1(1), p2) .* ...
	    	       sf_butter(intri(:, 2), p1(2), p2);
	    	extro = extro .* filt(:, [1 1 1]);
    	    else
    	    	error(['Invalid parameter for ', aj, '.']);
	    end
    	
    	
    	%%% anything else is a mystery
    	otherwise,
    	    error(['Unknown filter: ', aj, '.']);
    end
    j = j+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting function(s)

function butter = sf_butter(freq, bandwidth, order)
% SF_BUTTER	Compute Butterworth filter transfer function
%		butter = SF_BUTTER(freq, bandwidth, order)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin ~= 3) 
    help sf_butter; return;
elseif prod(size(bandwidth)) ~= 1
    error(['Bandwidth must be a scalar.']);
elseif prod(size(order)) ~= 1
    error(['Filter must be a scalar.']);
end

p = i * freq / bandwidth;

switch order
  case 1, 	butter = 1 ./ (1 + p)
  case 2,	butter = 1 ./ (1 + sqrt(2) * p + p.^2)
  case 3,	butter = 1 ./ ((1 + p) .* (1 + p + p.^2))
  case 4,	butter = 1 ./ ((1 + 0.765 * p + p.^2) .* (1 + 1.848 * p + p.^2))
  case 5,	butter = 1 ./ ((1 + p) .* (1 + 0.618 * p + p.^2) .* (1 + 1.62 * p + p.^2))
  otherwise,	error(['Unknown or unimplemented filter order ', num2str(order), '.'])
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
