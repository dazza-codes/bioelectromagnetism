function [ointr, er] = sf_opt_intr(method, limit, lval, intr, extr, border, edge, face)
% SF_OPT_INTR	Optimize intrinsic coordinates
%		[ointr, er] = SF_OPT_INTR(method, limit, lval, intr, extr, border, edge [, face])
%		use METHOD, LIMIT, LVAL to compute optimized points OINTR
%   	    	from given INTRs, EXTRs, BORDER, EDGEs, FACEs
%		ER = optional output for #pass/time/mflops/errd/errl/deld/dell
%		set global SF_OPT_INTR_POW to exponent (default == 1)
%	METHOD = { 'difc', 'dife', 'diff', 'fmins', 'none' }
%		appended 'v' displays info at each iteration
%		appended 'q' suppresses all info
%	    DIFC:  diffusion from edge connectivity
%	    DIFE:  diffusion from edge length
%	    DIFF:  diffusion from face area
%	    FMINS: MATLAB function
%	    NONE:  no optimization (base case) 
%	LIMIT = { 'deld', 'dell', 'errd', 'errl', 'flop', 'pass', 'time' }
%	    DELD: delta difference  error between successive iterations
%	    DELL: delta logarithmic error between successive iterations
%	    ERRD: difference  error (computed & normed by sf_opt_intr_mf)
%	    ERRL: logarithmic error (computed & normed by sf_opt_intr_mf)
%	    FLOP: number of Mflops
%	    PASS: number of passes
%	    TIME: elapsed time

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - according to mprof
%    85% flops and 70% time in sf_opt_intr_mf
%    10% flops and  5% time in ointr = diffmat*ointr
%
% - errors often increase early or occasionally in optimization
% - to get final #pass/time/mflops:	er(size(er,1),1:3)
% - approaches:
%   - iterative (self-averaging) diffusion
%     - could formulate with relaxation parameter, 
%       although tests suggest that optimal value is close to 1
%     - edge cnnctvty	- results in uniform edge length 
%			  (and hence uniform area, uniform density)
%     - edge length	- matches intrinsic and extrinsic densities
%     - face area	- matches intrinsic and extrinsic densities
%   - MATLAB optimization procedures (neither very successful so far)
%     - FMINS (use exponents 1,2,3?)

%%% THINGS TO DO
% - get FMINS working (exponents 1,2,3)
% ? extra output argument for diffusion matrix
% ? use PROPAGATION SPEED and size to determine # passes
% ? rescale points to full range at each pass instead of using border
%   - may have tried this without luck - corners drift in?
% ? attraction/repulsion spacing force (Szeliski etal)  
% ? require all inputs to be real (faster)
% ? increase diffusion efficiency - see Press etal pg 652ff
%   ? adjust relaxation parameter (currently commented)
%   ? compare to parameters in sf_def_surf...
% ? wrap edges/faces?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 7)
    help sf_opt_intr; return; 
elseif ~ischar(method)
    error(['METHOD must be a string.']);
elseif ~ischar(limit)
    error(['LIMIT must be a string.']);
elseif prod(size(lval)) ~= 1
    error(['LVAL must be a scalar.']);
elseif (size(intr, 2) ~= 2)
    error(['INTR must have 2 columns,  not ', num2str(size(intr, 2)), '.']);
elseif all(size(extr, 2) ~= [3 2])
    error(['EXTR must have 3 (or 2) columns,  not ', num2str(size(extr, 2)), '.']);
elseif size(intr, 1) ~= size(extr, 1)
    error(['INTR and EXTR must have the same number of rows.']);
elseif size(border, 2) ~= 1
    error(['BORDER must have 1 columns,  not ', num2str(size(border, 2)), '.']);
elseif size(border, 1) ~= size(intr, 1)
    error(['BORDER and INTR must have the same number of rows.']);
elseif size(edge, 2) ~= 2
    error(['EDGE must have 2 columns,  not ', num2str(size(edge, 2)), '.']);
end
if (nargin > 7) 
    if size(face, 2) ~= 3
	error(['FACE must have 3 columns, not ', num2str(size(face, 2)), '.']);
    end
end

%%% determine output verbosity
dval = 1;
if     (method(end) == 'v')	dval = 2;
elseif (method(end) == 'q')	dval = 0;
end

if (dval > 0)
    dformat = ' %3d %7.3f %7.3f %7.3f %7.2f %7.3f %7.3f %s';
    disp([' ']);
    disp(['pass    time   Mflop   errd    errl   /\ d    /\ l  |method limit  lval  ']);
    disp(['----+-------+-------+-------+-------+-------+-------|------+-----+-------']);
    label = sprintf(' %4s %4s %7.3f \n', upper(method), upper(limit), lval);
end

power  	= 1;

global SF_OPT_INTR_POW
if ~isempty(SF_OPT_INTR_POW) power = SF_OPT_INTR_POW; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% diffusion methods

if strcmp(method(1:3), 'dif')
    
    flop1 = flops;
    time1 = cputime;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% diffusion based on edge connectivity

    if strcmp(method(1:4), 'difc')
	
	% get row, col, val entries from each edge
	%   (row,col) are indices for (val) in sparse diffusion matrix
	row = reshape(edge(:, [1, 2]), prod(size(edge)), 1);
	col = reshape(edge(:, [2, 1]), prod(size(edge)), 1);
	val = ones(size(row, 1), 1);
	

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% diffusion based on extrinsic edge length
    %%%   closer points (smaller edge lengths) should have more weight
    elseif strcmp(method(1:4), 'dife')
	
	% get row, col, val entries from each edge
	%   (row,col) are indices for (val) in sparse diffusion matrix
	epts = reshape(extr(edge,:),[size(edge) size(extr,2)]);
	val  = 1 ./ sqrt(sum(diff(epts,1,2).^2,3));
	val = [val; val];
	row = reshape(edge(:, [1 2]), prod(size(edge)), 1);
	col = reshape(edge(:, [2 1]), prod(size(edge)), 1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% diffusion based on extrinsic face area
    %%%    closer points (smaller face areas) should have more weight
    elseif strcmp(method(1:4), 'diff')
    
	% get row, col, val entries from each face
	%   (row,col) are indices for (val) in sparse diffusion matrix
	fpts = reshape(extr(face(:,[1 2 3 1]),:),[size(face,1) 4 size(extr,2)]);
	flen = sqrt(sum(diff(fpts,1,2).^2,3));
	frad = sum(flen,2)/2;
	val = 1 ./ sqrt(prod([frad, frad(:,[1 1 1]) - flen],2));
	val = [ val, val, val, val, val, val ];
	row = reshape(face(:, [1, 1, 2, 2, 3, 3]), prod(size(face))*2, 1);
	col = reshape(face(:, [2, 3, 1, 3, 1, 2]), prod(size(face))*2, 1); 


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% invalid diffusion method
    else
	error(['Invalid METHOD ''', method, '''.']);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create sparse diffusion matrix 
    %   (diag values for borders, off-diag weights for non-borders)
    
    spndx = 1:length(border);
    diffmat = sparse(spndx,spndx,  border) + ...
    	   sparse(spndx,spndx, ~border) * ...
	     sparse(row,col, val, size(intr, 1),size(intr, 1));
  
    % normalize so that each row sum is 1 (for c = D * c)
    diffmat = sparse(spndx,spndx, 1 ./ sum(diffmat,2)) * diffmat;

    % relaxation parameter 
    % (results improve from 0 to 1, but may become unstable thereafter)
    % diffmat = (1-lambda) * sparse(spndx,spndx,1) + ...
    %	     (  lambda) * diffmat;
    
    flop_count = flops - flop1;
    time_count = cputime - time1;

    %%% temporary values
    pass	= 1;
    ointr	= intr;

    oerrd	= sf_opt_intr_mf(ointr, extr, edge, 'errd', power);
    oerrl	= sf_opt_intr_mf(ointr, extr, edge, 'errl', power);
    er = [ 0, time_count, flop_count/1000000, oerrd, oerrl, oerrd, oerrl ];
    if (dval == 2) fprintf(1, dformat, er, label); end
    
    %%% start iterative diffusion
    while (1)
	%%% perform iteration
	flop1 = flops;
	time1 = cputime;
	ointr = diffmat * ointr;
	flop_count = flop_count + flops - flop1;
	time_count = time_count + cputime - time1;

%    	sf_xplot(ointr, extr, face); pause
    	
	%%% compute error (not required)
	errd = sf_opt_intr_mf(ointr, extr, edge, 'errd', power);
	errl = sf_opt_intr_mf(ointr, extr, edge, 'errl', power);

	er   = [ er; pass, time_count, flop_count/1000000, ...
		 errd, errl, oerrd - errd, oerrl - errl];
	pass = pass + 1;
	if (dval == 2) fprintf(1, dformat, er(pass, :), label); end
	
	%%% abort if limit is reached (errors may increase initially)
	switch limit
	    case 'deld', if lval > abs(oerrd - errd),	break; end
	    case 'dell', if lval > abs(oerrl - errl),	break; end
	    case 'errd', if lval > errd,    	    	break; end
	    case 'errl', if lval > errl,    	    	break; end
	    case 'flop', if lval < flop_count/1000000,	break; end
	    case 'pass', if lval < pass,    	    	break; end
	    case 'time', if lval < time_count,	    	break; end
	    otherwise, error(['unknown limit: ',limit,'.']);
	end
	if (dval > 0) & (errd > oerrd) & (errl > oerrl)
	    disp(['WARNING: sf_opt_intr errors increasing at pass ', num2str(pass)]);
	end

	oerrd = errd; oerrl = errl;
    end
    if (dval == 1) fprintf(1, dformat, er(pass, :), label); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% null method

elseif strcmp(method(1:4), 'none') | strcmp(method(1:4), 'null')

    ointr = intr;

    errd = sf_opt_intr_mf(ointr, extr, edge, 'errd', power);
    errl = sf_opt_intr_mf(ointr, extr, edge, 'errl', power);

    er = [ 0, 0, 0, errd, errl, errd, errl ];
    if (dval > 0) fprintf(1, dformat, er, label); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB minimization methods

elseif strcmp(method(1:5), 'fmins')

    flop1 = flops;
    time1 = cputime;

    options = foptions;

    ointr = fmins('sf_opt_intr_mf', intr, options, [], ...
		  extr, edge, limit, power);

    flop_count = flops - flop1;
    time_count = cputime - time1;

    oerrd = sf_opt_intr_mf( intr, extr, edge, 'errd', power);
    oerrl = sf_opt_intr_mf( intr, extr, edge, 'errl', power);
     errd = sf_opt_intr_mf(ointr, extr, edge, 'errd', power);
     errl = sf_opt_intr_mf(ointr, extr, edge, 'errl', power);

    er = [ 0, time_count, flop_count/1000000, ...
	   errd, errl , oerrd - errd, oerrl - errl];
    if (dval > 0) fprintf(1, dformat, er, label); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% invalid method

else
    error(['Invalid METHOD ''', method, '''.']);
end

if (dval > 0) disp([' ']); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting function

function f = sf_opt_intr_mf(intr, extr, edge, method, power)
% SF_OPT_INTR_MF	merit function for optimizing intrinsic coordinates
%			f = SF_OPT_INTR_MF(intr, extr, edge, method, power)
%	METHOD = { 'errd', 'errl' }
%		ERRD: difference error
%		ERRL: logarithm  error

%%% NOTES
% - according to mprof:
%   55% - computing ilen/elen
%   30% - normalizing ilen/elen
%   15% - computing f
%
% - intrinsic and extrinsic edge lengths are normalized
%   so that the sum of each is 1.0 

%%% THINGS TO DO
% - USE LENGTHS, NOT SQUARED LENGTHS (not so important for diffusion)
% ? require all inputs to be real (faster)
% - add additional limits (and lvals)
% ? how to deal with 0-length intrinsic/extrinsic edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid and compatible arguments
if (nargin < 5)
    help sf_opt_intr_mf; return;
elseif all(size(intr, 2) ~= [2 3])
    error(['INTR must have 2 (or 3) columns, not ', num2str(size(intr, 2)), '.']);
elseif all(size(intr, 2) ~= [3 2])
    error(['EXTR must have 3 (or 2) columns, not ', num2str(size(extr, 2)), '.']);
elseif (size(intr, 1) ~= size(extr, 1))
    error(['INTR and EXTR must have the same length.']);
elseif (size(edge, 2) ~= 2)
    error(['EDGE must have 2 columns,  not ', num2str(size(edge, 2)), '.']);
elseif ~ischar(method)
    error(['METHOD must be a string.']);
elseif prod(size(power)) ~= 1
    error(['POWER must be a scalar.']);
end

%%% compute intrinsic and extrinsic lengths

ilen = sum((intr(edge(:, 1), :) - intr(edge(:, 2), :)).^2, 2);
elen = sum((extr(edge(:, 1), :) - extr(edge(:, 2), :)).^2, 2);

%%% normalize for comparison (add epsilon to avoid numeric problems below)
ilen = eps + ilen / sum(ilen);
elen = eps + elen / sum(elen);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% pick specific method

if     strcmp(method, 'errd')	f = sum(abs(    ilen -  elen) .^power);
elseif strcmp(method, 'errl')	f = sum(abs(log(ilen ./ elen)).^power);
else
    error(['Invalid METHOD ''', method, '''.']);
end

if f == 0
    disp(['WARNING: sf_opt_intr_mf returning 0 - invoking keyboard...']);
    keyboard;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
