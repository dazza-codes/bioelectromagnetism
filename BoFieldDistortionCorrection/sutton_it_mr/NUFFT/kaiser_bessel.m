 function [kb, alpha, kb_m] = kaiser_bessel(x, J, alpha, kb_m, K_N)
%function [kb, alpha, kb_m] = kaiser_bessel(x, J, alpha, kb_m)
%function [kb, alpha, kb_m] = kaiser_bessel(x, J, 'best', 0, K_N)
%
%	generalized Kaiser-Bessel function for x in support [-J/2,J/2]
%	shape parameter "alpha" (default 2.34 J)
%	order parameter "kb_m" (default 0)
%	see (A1) in lewitt:90:mdi, JOSA-A, Oct. 1990
%	in
%		x	[M,1]	arguments
%	out
%		kb	[M,1]	KB function values, if x is numbers
%				or string for kernel(k,J), if x is 'string'
%				or inline function, if x is 'inline'
%		alpha
%		kb_m
%
%	Copyright 2001-3-30	Jeff Fessler	The University of Michigan

%	if no arguments, make example plots
if nargin < 2
	help(mfilename)
	J = 8; alpha = 2.34 * J;
	x = linspace(-(J+1)/2, (J+1)/2, 1001)';
%	x = linspace(J/2-1/4, J/2+1/4, 1001)';

	mlist = [-4 0 2 7];
	leg = {};
	for ii=1:length(mlist)
		kb_m = mlist(ii);
		yy(:,ii) = kaiser_bessel(x, J, alpha, kb_m);
		func = kaiser_bessel('inline', 0, alpha, kb_m);
		yf = func(x, J);
		if any(yf ~= yy(:,ii)),
		[yf yy(:,ii)]
		error 'bug', end
		leg{ii} = sprintf('m=%d', kb_m);
	end
	yb = kaiser_bessel(x, J, 'best', [], 2);
	plot(	x, yy(:,1), 'c-', x, yy(:,2), 'y-', ...
		x, yy(:,3), 'm-', x, yy(:,4), 'g-', x, yb, 'r--')
	leg{end+1} = 'best';
	axis tight, legend(leg)
%	axisy(0, 0.01)	% to see endpoints
	xlabel \kappa, ylabel F(\kappa)
	title(sprintf('KB functions, J=%g \\alpha=%g', J, alpha) )
return
end

if ~isvar('J'), J = 6; end
if ~isvar('alpha') | isempty('alpha'), alpha = 2.34 * J; end
if ~isvar('kb_m') | isempty('kb_m'), kb_m = 0; end

if ischar(alpha)
	if streq(alpha, 'best')
		if K_N ~= 2, error 'only K/N=2 done', end
		s = 'private/kaiser,m=0';
		kb_m = 0;	% hardwired, because it was nearly the best!
                s = load(s);
		ii = find(J == s.Jlist);
		if isempty(ii)
			ii = imin(abs(J - s.Jlist));
			warning(sprintf('J=%d not found, using %d', ...
				J, s.Jlist(ii)))
		end
		alpha = J * s.abest.zn(ii);
	else
		error 'unknown alpha mode'
	end
end

if ischar(x)
	if ischar(alpha)
		if ~isvar('K_N'), error 'K_N required', end
		kb = 'kaiser_bessel(k, J, ''%s'', [], %g)';
		kb = sprintf(kb, alpha, K_N);
	else
		kernel_string = 'kaiser_bessel(k, J, %g, %g)';
		kb = sprintf(kernel_string, alpha, kb_m);
	end
	if streq(x, 'inline')
		kb = inline(kb, 'k', 'J');
	elseif ~streq(x, 'string')
		error '1st argument must be "inline" or "string"'
	end
return
end

ii = abs(x) < J/2;
f = sqrt(1 - (x(ii)/(J/2)).^2);
denom = besseli(kb_m,alpha);
if ~denom
	Sprintf('m=%g alpha=%g', kb_m, alpha)
end
kb = zeros(size(x));
kb(ii) = f.^kb_m .* besseli(kb_m, alpha*f) / denom;
kb = Real(kb);
