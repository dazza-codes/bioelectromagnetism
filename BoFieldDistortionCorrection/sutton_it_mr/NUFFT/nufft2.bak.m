 function X = nufft2(x, st)
%function X = nufft2(x, st)
%	Compute 2D NUFFT of image x
%	in:	x	[N1,N2,L]	input images
%		st			structure precomputed by nufft2_init()
%	out:	X	[M,L]		output spectra
%
%	Copyright 2001-9-17	Jeff Fessler	The University of Michigan

%
%	if no arguments, then run a simple test
%
if nargin < 1
	help(mfilename)

	J1 = 5; J2 = 6;
	N1 = 60; N2 = 75;
	K1 = 2*N1;	gam1 = 2*pi/K1;
	K2 = 2*N2;	gam2 = 2*pi/K2;
	n_shift = [0 0];
	Sprintf('err alf1 %g best %g', ...
		max(col(nufft2_err('all', N1, N2, J1, J2, K1, K2, 1))), ...
		max(col(nufft2_err('all', N1, N2, J1, J2, K1, K2, 'best'))) )

	x = [[1:N1]'*ones(1,3), ones(N1,N2-3)]; % test signal
%	x = randn(N1,N2);
%	x = zeros(N1,N2); x(1,1) = 1;
	if 0	% test with uniform frequency locations
		o1 = 2 * pi * [0:(N1-1)]' / N1;
		o2 = 2 * pi * [0:(N2-1)]' / N2;
		[o1, o2] = ndgrid(o1, o2);
		Xf = fft2(x);
	else	% nonuniform frequencies
		o1 = [0 7.2 2.6 3.3];
		o2 = [0 4.2 -1 5.5];
		[o1, o2] = ndgrid(linspace(0,gam1,11), linspace(0,gam2,13));
		om = [o1(:) o2(:)];
		Xd = dtft2(x, om, n_shift);
	end
	s0 = nufft2_init(om, N1, N2, J1, J2, K1, K2, n_shift, 0);
	s1 = nufft2_init(om, N1, N2, J1, J2, K1, K2, n_shift, 1);
	Sprintf('loop max %% difference = %g', max_percent_diff(s0.p,s1.p))
	Xn = nufft2(x, s1);
	Sprintf('alf1 max %% difference = %g', max_percent_diff(Xd,Xn))

	sb = nufft2_init(om, N1, N2, J1, J2, K1, K2, n_shift, 0, 'best');
	Xb = nufft2(x, sb);
	Sprintf('best max %% difference = %g', max_percent_diff(Xd,Xb))
return
end

N1 = st.N1;
N2 = st.N2;
K1 = st.K1;
K2 = st.K2;

dims = size(x);
if dims(1) ~= N1 | dims(2) ~= N2, error size, end

if 0 %round(K1/N1) == K1/N1 & round(K2/N2) == K2/N2
	%disp('note in nufft2: could save flops via smarter padded FFT')
        disp('saving flops via smarter padded FFT')
end

if ndims(x) > 2
	error 'not done'
else
	x = x .* st.sn;		% apply scaling factors
end

if 0 %round(K1/N1) == K1/N1 & round(K2/N2) == K2/N2
     % Not quicker in cases tried   03/18/02
        X1 = K1/N1;
        X2 = K2/N2;
        ind1 = [0:N1-1]';
	ind2 = [0:N2-1];
        for ii = 1:X1
            for jj = 1:X2
        	Xk(ii:X1:K1,jj:X2:K2) = fft2(x.*(exp(-i*2*pi/K1*(ii-1)*(ind1))*exp(-i*2*pi/K2*(jj-1)*(ind2))), N1, N2);
	    end
	end
else
	Xk = fft2(x, K1, K2);	% [K1,K2,L] oversampled FFT, padded at end
end

if ndims(x) > 2
	Xk = reshape(Xk, [K1*K2 dims(3:end)]);	% [K1*K2,L]
else
	Xk = Xk(:);
end

%
%	interpolate using precomputed sparse matrix
%
X = st.p * Xk;
