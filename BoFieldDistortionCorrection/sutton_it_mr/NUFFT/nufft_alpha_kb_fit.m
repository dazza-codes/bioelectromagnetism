 function [alphas, beta] = nufft_alpha_kb_fit(N, J, K, L, beta, chat)
%
%	return the alpha and beta corresponding to LS fit of L components
%	to optimized Kaiser-Bessel scaling factors (m=0, alpha=2.34J).
%	This is the best method I know currently for choosing alpha!
%
%	Copyright 2002-7-16	Jeff Fessler	The University of Michigan

if nargin < 3, help(mfilename), error(mfilename), end
if ~isvar('L') | isempty(L),		L = 13; end
if ~isvar('beta') | isempty(beta),	beta = 1; end
if ~isvar('chat') | isempty(chat),	chat = 0; end

kb_alf = 2.34 * J;	% KB shape parameter
kb_m = 0;		% KB order

[tmp, sn_kaiser] = nufft1_error(0, N, J, K, 'kaiser', 'ft');
sn_kaiser = Real(sn_kaiser);

%
% use regression to match NUFFT with BEST kaiser scaling's
%
gam = 2*pi/K;
nlist = [0:(N-1)]' - (N-1)/2;
X = cos(beta*gam*nlist*[0:L]);
coef = regress(sn_kaiser, X)';
alphas = [Real(coef(1)) coef(2:end)/2];

if chat
	Sprintf('cond # for LS fit to KB scale factors: %g', cond(X))
end
