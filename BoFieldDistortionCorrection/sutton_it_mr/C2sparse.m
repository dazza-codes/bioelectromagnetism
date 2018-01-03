 function [C, wjk] = C2sparse(type, kappa, nbrs, chat)
%function [C, wjk] = C2sparse(type, kappa, nbrs, chat)
%
%	create 2D penalty matrix C with 2nd-order pixel neighborhood
%	R = beta * C'*D(wjk)*C is the penalty Hessian (in quadratic case)
%	in:
%		type	string		choices:
%					'leak' - allow penalty across mask edge
%						(consistent with aspire)
%					'tight' - no penalty across mask edges
%		kappa	[n1,n2]		if logical, then just support (aka mask)
%					if double, then as in fessler:96:srp
%		nbrs	scalar		4 or 8 neighbors
%	out:
%		C [n1*n2*nbrs/2,n1*n2]	most rows have a 1 and a -1
%		wjk [n1*n2*nbrs/2,1]	penalty weights:
%					kappa_j * kappa_k * [1 or 1/sqrt(2)]
%
%	Caller may want to do:	C = C(:,mask(:));
%
%	Copyright 2002-1-30	Jeff Fessler	The University of Michigan

%
%	default is to give help and run a test example
%
if nargin < 3 | ~ischar(type)
	help(mfilename)
	n1 = 5; n2 = 4;
	kappa = ones(n1,n2);
	hood = 8;
	[C wjk] = C2sparse('leak', kappa, hood, 1);
return
end


if ~isvar('chat'), chat = 0; end

if streq(type, 'leak')
	leak = logical(1);
elseif streq(type, 'tight')
	leak = logical(0);
else
	error 'bad type'
end

if nbrs == 8
	offset1 = [-1 +0 -1 +1];
	offset2 = [+0 -1 -1 -1];
elseif nbrs == 4
	offset1 = [-1 +0];
	offset2 = [+0 -1];
else
	error 'bad nbrs'
end

[n1 n2 n3] = size(kappa);
if n3 > 1, error '3d not done', end


%
%	build matrix, stacking up a block for each neighbor
%
j1 = 1:n1;
j2 = 1:n2;
[j1 j2] = ndgrid(j1, j2);
jj = sub2ind([n1 n2], j1, j2);

C = [];
wjk = zeros(n1*n2, length(offset1));

%
%	loop over all neighboring pixels
%
for io=1:length(offset1)
	o1 = offset1(io);
	o2 = offset2(io);

	k1 = j1 + o1;	% k: index of neighbor
	k2 = j2 + o2;

	%
	% only make '-1' entries for k's that are within rectangle AND mask
	%
	krect = k1 >= 1 & k1 <= n1 & k2 >= 1 & k2 <= n2;
	kk = ones(size(k1));	% dummy '1' value for non-rect cases
	kk(krect) = sub2ind([n1 n2], k1(krect), k2(krect));

	mask = kappa(:) > 0;
	kmask = zeros(size(mask));
	kmask(krect(:)) = mask(kk(krect));

	if leak		% put 1 or -1 if either j or k is within mask
		jvalue = mask;
		kvalue = kmask;
	else		% both j and k must be within mask
		jvalue = mask & kmask;
		kvalue = jvalue;
	end
	Ctmp	= sparse(1:n1*n2, jj(:), jvalue, n1*n2, n1*n2) ...
		- sparse(1:n1*n2, kk(:), kvalue, n1*n2, n1*n2);

	wjk(:,io) = (jvalue | kvalue) / sqrt(o1^2 + o2^2);

	%
	%	include kappa effect, but handle "leak" kappas carefully
	%	cf: wjk(:,io) *= kappa(jj) .* kappa(kk)
	%
	kappa2 = zeros(size(mask));
	jk = mask & kmask;
	kappa2(jk) = kappa(jk) .* kappa(kk(jk));
	jk = mask & ~kmask;
	kappa2(jk) = kappa(jk).^2;
	jk = ~mask & kmask;
	kappa2(jk) = kappa(kk(jk)).^2;
	wjk(:,io) = wjk(:,io) .* kappa2(:);

	C = [C; Ctmp];
end

wjk = wjk(:);

%
%	optional display
%
if chat
	clf
	subplot(131), spy(C), title('C')
	if ncol(C) < 500
		im(132, C'*C, 'C''C')
	else
		im(132, embed(diag(C'*C), kappa>0), 'Diag(C''C)'), cbar
	end
	tmp = reshape(wjk, [n1 n2 nbrs/2]);
	im(133, tmp, 'wjk'), cbar
end
