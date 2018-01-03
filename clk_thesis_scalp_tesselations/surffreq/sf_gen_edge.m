function edge = sf_gen_edge(intr, max)
% SF_GEN_EDGE	Generate set of non-intersecting edges for intrinsic points
%		edge = SF_GEN_EDGE(intr [, max])
%		ignores potential edges with length > MAX
%		(work done by C MEX-file)
% adj matrix:	sparse([ed(:,1); ed(:,2)], [ed(:,2); ed(:,1)], 1)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

