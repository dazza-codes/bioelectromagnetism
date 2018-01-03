function extr = sf_gen_extr(intr, shape, par, range)
% SF_GEN_EXTR	Generate set of 3D extrinsic locations
%		extr = SF_GEN_EXTR(intr, shape [, par [, range]])
%		  SHAPE = {'cos_u', 'cos_v', 'cos_uv', 'cylinder', 
%			'impulse', 'sheet', 'sin_u', 'sin_v', 'sin_uv', 
%			'sinc_u', 'sinc_v', 'sinc_uv', 'sphere', 'torus'}
%		using PAR as a control value
%		maps each INTRinsic point to (-RANGE/2 ... RANGE/2)^2
%		defaults:  RANGE=1.0  PAR=1.0

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO
% ? other shapes - 4gon, 6gon, 8gon
% ? more control over surface orientation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_gen_extr; return; end
if (nargin < 3) par = 1.0; end
if (nargin < 4) range = 1.0; end

extr = zeros(size(intr, 1), 3);

if size(intr, 2) ~= 2
    error(['Invalid number of columns (', num2str(size(intr, 2)), ') in INTR.']);
elseif ~ischar(shape)
    error(['SHAPE must be a string.']);
elseif prod(size(range)) ~= 1
    error(['RANGE must be a scalar.']);
elseif prod(size(par)) ~= 1
    error(['PAR must be a scalar.']);
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch shape

    case 'cos_u',
    	extr = [ intr, 0.5 * cos(par*2*pi*intr(:,1)) ];
    case 'cos_v',
	extr = [ intr, 0.5 * cos(par*2*pi*intr(:,2)) ];
    case 'cos_uv',
	extr = [ intr, 0.5 * cos(par*2*pi*intr(:,1)) .* cos(par*2*pi*intr(:,2)) ];

    case 'cylinder',
	extr = [ 0.5 * sin(2*pi*intr(:,1)), ...
	    	 0.5 * cos(2*pi*intr(:,1)), ...
	    	 0.5 * par *  2*intr(:,2)];

    case 'impulse',
	extr = [ intr, zeros(size(intr, 1), 1) ];
	extr((extr(:,1).^2 + extr(:,2).^2 < 0.01 * par),3) = 0.5;
	
    case 'sheet',
	extr = [ intr, zeros(size(intr, 1), 1) ];

    case 'sin_u',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,1)) ];
    case 'sin_v',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,2)) ];
    case 'sin_uv',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,1)) .*  sin(par*2*pi*intr(:,2)) ];

    case 'sinc_u',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,1)) ./ (par*2*pi*intr(:,1)) ];
    case 'sinc_v',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,2)) ./ (par*2*pi*intr(:,2)) ];
    case  'sinc_uv',
	extr = [ intr, 0.5 * sin(par*2*pi*intr(:,1)) .* sin(par*2*pi*intr(:,2)) ./ ...
			((par*2*pi*intr(:,1)) .*    (par*2*pi*intr(:,2))) ];

    case 'sphere',
	extr = [ 0.5 * sin(2*pi*intr(:,1)) .* cos(pi*intr(:,2)), ...
		 0.5 * cos(2*pi*intr(:,1)) .* cos(pi*intr(:,2)), ...
		 0.5 * par *                  sin(pi*intr(:,2))];

    case 'torus',
	extr = [ 0.5 * sin(2*pi*intr(:,1)) .* (1.0 + 0.5 * cos(2*pi*intr(:,2))), ...
		 0.5 * cos(2*pi*intr(:,1)) .* (1.0 + 0.5 * cos(2*pi*intr(:,2))), ...
		 0.5 * par *                               sin(2*pi*intr(:,2))];

    otherwise,
	error(['Unknown SHAPE ''', shape, '''.']);
end

extr = range*extr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
