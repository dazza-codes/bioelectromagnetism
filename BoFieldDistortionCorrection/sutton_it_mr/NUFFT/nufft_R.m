 function R = nufft_R(N, J, K, tol, alpha, gamma_scale)
%function R = nufft_R(N, J, K, tol, alpha, gamma_scale)
%	Precompute the matrix R = [C' S S' C]\inv used in NUFFT.
%	This can be precomputed, being independent of frequency location.
%	in:
%		N		# signal length
%		J		# of neighbors
%		K		# FFT length
%		tol		tolerance for smallest eigenvalue
%		alpha	[L+1]	Fourier coefficient vector for scaling
%		gamma_scale	scale gamma=2*pi/K by this for Fourier series
%	out:
%		R	[J,J]	precomputed matrix
%
%	Copyright 2000-1-9	Jeff Fessler	The University of Michigan
%
%	if no arguments, run a test
%
if nargin < 3
	help(mfilename)

	N = 1; K = 2*N;
	alpha = [1 0 0]
	for J=1:8
		R = nufft_R(N, J, K, [], alpha);
		Sprintf('J=%d K/N=%d cond=%g', J, K/N, cond(R))
	end
	clear R, return
end

if ~isvar('tol') | isempty(tol)
	tol = 1e-7;
end

if ~isvar('gamma_scale') | isempty(gamma_scale)
	gamma_scale = 1/2;
end

if N > K, error 'N > K', end


%
%	default with unity scaling factors
%
if ~isvar('alpha') | isempty(alpha)

	%
	%	compute C'SS'C = C'C
	%
	[j1, j2] = ndgrid(1:J, 1:J);
	cssc = nufft_diric(j2 - j1, N, K);


%
%	Fourier-series based scaling factors
%
else
	if ~isreal(alpha(1)), error 'need real alpha_0', end
	L = length(alpha) - 1;	% L
	cssc = zeros(J,J);
	[j1, j2] = ndgrid(1:J, 1:J);
	for l1 = -L:L
		for l2 = -L:L
			alf1 = alpha(abs(l1)+1);
			if l1 < 0, alf1 = conj(alf1); end
			alf2 = alpha(abs(l2)+1);
			if l2 < 0, alf2 = conj(alf2); end

			tmp = j2 - j1 + gamma_scale * (l1 - l2);
			tmp = nufft_diric(tmp, N, K);
			cssc = cssc + alf1 * conj(alf2) * tmp;
%		Sprintf('%d %d %s %s', l1, l2, num2str(alf1), num2str(alf2))
		end
	end
end


%
%	Inverse, or, pseudo-inverse
%

%smin = svds(cssc,1,0);
smin = min(svd(cssc));
if smin < tol	% smallest singular value
	warning(sprintf('Poor conditioning %g => pinverse', smin))
	R = pinv(cssc, tol/10);
else
	R = inv(cssc);
end
