function [S, er] = sf_xform(s, e, E, method, limit, lval, p1)
% SF_XFORM	Transform between spatial and frequency domains
%		[S, er] = SF_XFORM(s, e, E, method [, limit, lval [, p1]])
%		computes output coordinates S
%		from input coordinates s
%		using coefficients e (forward) and E (inverse)
%		ER = optional output for #pass/time/Mflops/e1/e2/ef/d1/d2/df
%	METHOD = { 'cor', 'dir', 'for', 'inv', 'tra' }
%	    appended 'v' displays info at each iteration
%	    appended 'q' suppresses all info
%	    COR: iterative correction (P1 = L)
%		S0 = Les, M = (I - LeE)		Sn+1 = S0 + MSn
%	    DIR: direct solution 		(ignores LIMIT)
%		S = es
%	    DIV: matrix division using MATLAB 	(ignores LIMIT)
%		S = E/s
%	    FOR: iterative forward (P1 = L)
%		f0 = Le, M0 = I, M = (I - LeE)	Mn+1 = M * Mn, fn+1 = fn + Mn+1
%	    INV: matrix inverse using MATLAB 	(ignores LIMIT)
%		S = inv(E)s
%	    TRA: iterative transformation (P1 = L)
%		S0 = Les	    		Sn+1 = Sn + Le( s - ESn)
%	LIMIT = { 'del1', 'del2', 'delf', 'err1', 'err2', 'errf', 
%                 'flop', 'pass', 'time' }
%	    DEL1: delta 1-norm error between successive iterations
%	    DEL2: delta 2-norm error between successive iterations
%	    DELF: delta F-norm error between successive iterations
%	    ERR1: 1-norm of error (divided by surface norm)
%	    ERR2: 2-norm of error (divided by surface norm)
%	    ERRF: F-norm of error (divided by surface norm)
%	    FLOP: number of Mflops
%	    NULL: no limit
%	    PASS: number of passes
%	    TIME: elapsed time

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - to get final #pass/time/mflops:	er(end,1:3)
% - associativity has been optimized for speed
% - need to change SPEYE to EYE to use MATLAB compiler

%%% THINGS TO DO
% - add spatial area for energy/power norms
% ? use default values for P1, limit, lval
% ? make limit and lval into arrays for multiple limits (how to break?)
% ? add freq area and redefine e and E?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here


%%% check for valid and compatible arguments
if (nargin < 4) 
    help sf_xform; return;
elseif size(s, 2) ~= 3
    error(['s must have 3 columns, not ', num2str(size(s,2)), '.']);
elseif size(e, 2) ~= size(s, 1)
    error(['e and s have incompatible sizes (', ...
           num2str(size(e, 2)), ' ~= ', num2str(size(s, 1)), ').']);
elseif size(E, 2) ~= size(e, 1)
    error(['E and e have incompatible sizes (', ...
	   num2str(size(E, 2)), ' ~= ', num2str(size(e, 1)), ').']);
elseif ~ischar(method) error(['METHOD must be a string.']);
elseif (nargin < 5)
    limit = 'null'; lval = 0;
elseif ~ischar(limit)  
    error(['LIMIT must be a string.']);
end

%%% determine output verbosity
dval = 1;
if size(method, 2) == 4
    if     strcmp(method(4), 'v')	dval = 2;
    elseif strcmp(method(4), 'q')	dval = 0;
    end
end

%%% temporary values
pass	= 1;
oerr1	= Inf; 
oerr2	= Inf; 
oerrf	= Inf;

snorm1	= norm(s, 1);
snorm2	= norm(s, 2);
snormf	= norm(s, 'fro');

if (dval > 0)
    dformat = ' %3d %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %s';
    disp([' ']);
    disp(['pass    time   Mflop   err1    err2    errf    /\ 1    /\ 2    /\ f |method limit  lval  ']);
    disp(['----+-------+-------+-------+-------+-------+-------+-------+-------|------+-----+-------']);
    label = sprintf(' %4s %4s %7.3f \n', upper(method), upper(limit), lval);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% direct solution (first approximation)
%%% matrix division
%%% matrix inverse

if strcmp(method(1:3), 'dir') | strcmp(method(1:3), 'div') | strcmp(method(1:3), 'inv')
    %%% generate estimate
    flop1 = flops;
    time1 = cputime;
    if     strcmp(method(1:3), 'dir')   S = e * s;
    elseif strcmp(method(1:3), 'div')	S = E \ s;
    elseif strcmp(method(1:3), 'inv')	S = inv(E) * s;
    end
    flop_count = flops - flop1;
    time_count = cputime - time1;

    %%% compute error (not required)
    err = s - E*S;
    err1 = norm(err, 1)		/ snorm1;
    err2 = norm(err, 2)		/ snorm2;
    errf = norm(err, 'fro')	/ snormf;

    er   = [ 0, time_count, flop_count/1000000, ...
	     err1, err2, errf, err1, err2, errf ];
    if (dval > 0) fprintf(1, dformat, er, label); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% iterative correction

elseif strcmp(method(1:3), 'cor')

    %%% check for valid and compatible arguments
    if (nargin ~= 7)
	error(['Invalid number of arguments (', num2str(nargin), ...
	       ') for METHOD ''cor''.']);
    end 
    
    flop1 = flops;
    time1 = cputime;
    S0 = e * (s * p1);
    S = S0;
    %%% WARNING - this takes a long time for full e and E
    m = speye(size(e, 1), size(E, 2)) - e * (E * p1);
    flop_count = flops - flop1;
    time_count = cputime - time1;

    er = [ 0, time_count, flop_count/1000000, ...
	   NaN, NaN, NaN, NaN, NaN, NaN ];
    if (dval == 2) fprintf(1, dformat, er, label); end
    
    %%% start iterative transform
    while (1)
	%%% generate estimate
	flop1 = flops;
	time1 = cputime;
	S = S0 + m * S;
	flop_count = flop_count + flops - flop1;
	time_count = time_count + cputime - time1;
	
	%%% compute error (not required)
	err = s - E*S;
	err1 = norm(err, 1)	/ snorm1;
	err2 = norm(err, 2)	/ snorm2;
	errf = norm(err, 'fro')	/ snormf;

	er   = [ er; pass, time_count, flop_count/1000000, ...
		     err1, err2, errf, ...
		     oerr1-err1, oerr2-err2, oerrf-errf ];
	pass = pass + 1;
	if (dval == 2) fprintf(1, dformat, er(pass, :), label); end

	%%% abort if limit is reached or all errors are increasing
	switch limit
	    case 'del1', if lval > abs(oerr1 - err1),	break; end
	    case 'del2', if lval > abs(oerr2 - err2),	break; end
	    case 'delf', if lval > abs(oerrf - errf),	break; end
	    case 'err1', if lval > err1,    	    	break; end
	    case 'err2', if lval > err2,    	    	break; end
	    case 'errf', if lval > errf,    	    	break; end
	    case 'flop', if lval < flop_count/1000000,	break; end
	    case 'pass', if lval < pass,    	    	break; end
	    case 'time', if lval < time_count,	    	break; end
	    otherwise, error(['unknown limit: ', limit, '.']);
	end
	if (err1 > oerr1) & (err2 > oerr2) & (errf > oerrf) 
	    disp(['WARNING: sf_xform errors increasing at pass ', num2str(pass)]);
	    break;
	end

	oerr1 = err1; oerr2 = err2; oerrf = errf;
    end
    if (dval == 1) fprintf(1, dformat, er(pass, :), label); end

	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% iterative forward

elseif strcmp(method(1:3), 'for')

    %%% check for valid and compatible arguments
    if (nargin ~= 7)
	error(['Invalid number of arguments (', num2str(nargin), ...
	       ') for METHOD ''for''.']);
    end 
    
    flop1 = flops;
    time1 = cputime;
    %%% LOOK OUT - this takes a long time for full e and E
    m = speye(size(e, 1), size(E, 2)) - e * (E * p1);
    mi = p1 * e;
    fwd = mi;
    % this line included so that flop_count and time_count are correct
    S = fwd * s;
    
    flop_count = flops - flop1;
    time_count = cputime - time1;
    
    er = [ 0, time_count, flop_count/1000000, ...
	   NaN, NaN, NaN, NaN, NaN, NaN ];
    if (dval == 2) fprintf(1, dformat, er, label); end

    %%% start iterative transform
    while (1)
	%%% generate estimate
	flop1 = flops;
	time1 = cputime;
   	%%% LOOK OUT - this takes FOREVER for full e and E
	mi = m * mi;
	fwd = fwd + mi;
	flop_count = flop_count + flops - flop1;
	time_count = time_count + cputime - time1;
	
	%%% compute error one way (not required)
	err = s - E*fwd*s;
	%%% compute error another way (not required)
%	err = E * fwd;
%	err = speye(size(err)) - err;

	err1 = norm(err, 1)	/ snorm1;
	err2 = norm(err, 2)	/ snorm2;
	errf = norm(err, 'fro')	/ snormf;

	er   = [ er; pass, time_count, flop_count/1000000, ...
		     err1, err2, errf, ...
		     oerr1-err1, oerr2-err2, oerrf-errf ];
	pass = pass + 1;
	if (dval == 2) fprintf(1, dformat, er(pass, :), label); end

	%%% abort if limit is reached or all errors are increasing
	switch limit
	    case 'del1', if lval > abs(oerr1 - err1),	break; end
	    case 'del2', if lval > abs(oerr2 - err2),	break; end
	    case 'delf', if lval > abs(oerrf - errf),	break; end
	    case 'err1', if lval > err1,    	    	break; end
	    case 'err2', if lval > err2,    	    	break; end
	    case 'errf', if lval > errf,    	    	break; end
	    case 'flop', if lval < flop_count/1000000,	break; end
	    case 'pass', if lval < pass,    	    	break; end
	    case 'time', if lval < time_count,	    	break; end
	    otherwise, error(['unknown limit: ', limit, '.']);
	end
	if (err1 > oerr1) & (err2 > oerr2) & (errf > oerrf) 
	    disp(['WARNING: sf_xform errors increasing at pass ', num2str(pass)]);
	    break;
	end

	oerr1 = err1; oerr2 = err2; oerrf = errf;
    end
    if (dval == 1) fprintf(1, dformat, er(pass, :), label); end
    S = fwd * s;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% iterative transformation

elseif strcmp(method(1:3), 'tra')

    %%% check for valid and compatible arguments
    if (nargin ~= 7)
	error(['Invalid number of arguments (', num2str(nargin), ...
	       ') for METHOD ''tra''.']);
    end 
    
    flop1 = flops;
    time1 = cputime;
% these represent a better "base case"
%    S	= zeros(size(e,1),size(err,2));
%    err	= s;
% these match other iterative methods
    S = e * (s * p1);
    err = s - E*S;
    flop_count = flops - flop1;
    time_count = cputime - time1;
    
    er = [ 0, time_count, flop_count/1000000, ...
	   NaN, NaN, NaN, NaN, NaN, NaN ];
    if (dval == 2) fprintf(1, dformat, er(1, :), label); end

    %%% start iterative transform
    while (1)
	%%% generate estimate using error
	flop1 = flops;
	time1 = cputime;
	S = S + e * (err * p1);

	%%% compute error (required)
	err = s - E*S;
	flop_count = flop_count + flops - flop1;
	time_count = time_count + cputime - time1;
	err1 = norm(err, 1)	/ snorm1;
	err2 = norm(err, 2)	/ snorm2;
	errf = norm(err, 'fro')	/ snormf;

	er   = [ er; pass, time_count, flop_count/1000000, ...
		     err1, err2, errf, ...
		     oerr1-err1, oerr2-err2, oerrf-errf ];
	pass = pass + 1;
	if (dval == 2) fprintf(1, dformat, er(pass, :), label); end

	%%% abort if limit is reached or all errors are increasing
	switch limit
	    case 'del1', if lval > abs(oerr1 - err1),	break; end
	    case 'del2', if lval > abs(oerr2 - err2),	break; end
	    case 'delf', if lval > abs(oerrf - errf),	break; end
	    case 'err1', if lval > err1,    	    	break; end
	    case 'err2', if lval > err2,    	    	break; end
	    case 'errf', if lval > errf,    	    	break; end
	    case 'flop', if lval < flop_count/1000000,	break; end
	    case 'pass', if lval < pass,    	    	break; end
	    case 'time', if lval < time_count,	    	break; end
	    otherwise, error(['unknown limit: ', limit, '.']);
	end
	if (err1 > oerr1) & (err2 > oerr2) & (errf > oerrf) 
	    disp(['WARNING: sf_xform errors increasing at pass ', num2str(pass)]);
	    break;
	end

	oerr1 = err1; oerr2 = err2; oerrf = errf;
    end
    if (dval == 1) fprintf(1, dformat, er(pass, :), label); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    error(['Unknown or unimplemented method (', method, ').']);
end

if (dval > 0) disp([' ']); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks



