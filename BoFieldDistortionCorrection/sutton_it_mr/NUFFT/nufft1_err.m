 function err = nufft1_err(om, N1, J1, K1, alpha, beta)
%function err = nufft1_err(om, N1, J1, K1, alpha, beta)
%	Compute worst-case error for each input frequency for 1D NUFFT.
%	in:
%		om	[M,1]	digital frequency omega in radians
%		N1		signal length
%		J1		# of neighbors used per frequency location
%		K1		FFT size (should be > N1)
%		alpha	[L,1]	Fourier series coefficients of scaling factors
%		beta		scale gamma=2pi/K by this in Fourier series
%				typically is K/N (me) or 0.5 (Liu)
%	out:
%		err	[M,1]	worst-case error over unit-norm signals
%
%	Copyright 2001-12-7	Jeff Fessler	The University of Michigan

%	if no arguments, give an example
if nargin < 4
	help(mfilename)
	N = 1; K = 2*N; J = 6; gam = 2*pi/K;
	om = gam * linspace(0,1,101);
	plot(om/gam, nufft1_err(om, N, J, K))
	xlabel '\omega / \gamma', ylabel 'E_{max}(\omega)'
	return
end

if ~isvar('alpha') | isempty(alpha)
	alpha = [1];	% default Fourier series coefficients of scaling factors
end
if ~isvar('beta') | isempty(beta)
	beta = 0.5;	% default is Liu version for now
end

%
%	see if 'best' alpha is desired
%
if ischar(alpha)
	if ~strcmp(alpha, 'best'), error 'unknown alpha argument', end
	[alpha, beta, ok] = nufft_best_alpha(J1, 0, K1/N1);
	if ~ok, error 'unknown J,K/N', end
end

tol = 0;
R1 = nufft_R(N1, J1, K1, tol, alpha, beta);	% [J,J]
r1 = nufft_r(om, N1, J1, K1, alpha, beta);	% [J,M]

%
%	worst-case error at each frequency
%
Rr1 = R1 * r1;			% [J,M]
err = sum(conj(r1) .* Rr1).';	% [M,1]
err = min(real(err), 1);
err = sqrt(1 - err);
