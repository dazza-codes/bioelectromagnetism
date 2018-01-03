function S = sf_nspec(c,C,a,A,s,xform_args)
% SF_NSPEC  Compute normalized spectrum (used in test_def*)
%   	    S = sf_nspec(c,C,a,A,s,xform_args)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

if (nargin < 6)
    help sf_nspec; return;
end

[map,fmap] = sf_gen_map(c,C,'ft');
 map =  map * sparse(1:length(a),1:length(a), a);
fmap = fmap * sparse(1:length(A),1:length(A), A);
S = sf_xform(sf_gxform(s,'nrm',[]), map,fmap, xform_args{:});

