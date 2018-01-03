function [intr, extr, brdr, edge, face, nrby, err] = sf_def_surf(intr, extr, brdr, edge, face, target, dterms)
% SF_DEF_SURF	Deform surface (extrinsic coordinates) to match target
%   	    	[intr, extr, brdr, edge, face, nrby, err] = ...
%   	    	    SF_DEF_SURF(intr, extr, brdr, edge, face, target, dterms)
%   	    	(surfaces with INTR/EXTR coordinates, BRDRs, EDGEs, FACEs)
%   	    	NRBY is logical vector of points that reached target
%   	    	TARGET is 2D/3D array of 3D data points on surface
%   	    	DTERMs in deformation equation (cell array):
%
%   	batch	- run in batch mode (no user interaction)
%
%   	pass	- specify number of iterations (dummy)	(# iterations)
%   	point	- specify number of points (dummy)  	(# points)
%   	passpoint - specify pass point product (dummy)	(product)
%
%   	ssize	- step size for other terms 	    	(method ,[ max step])
%   	    	    ('fixed', 'gmean', 'gmin', 'lmean')
%
%   	datas	- move toward closest scattered point	(step size, search dist, max num)
%   	datav  	- move toward nearby volumetric point	(step size, search dist)
%
%   	  dir* only affect subset of points: ('all','dstnt','nrby','none')
%   	dirr	- move in random direction		(step size, subset)
%   	dirn	- move along face surface normals	(step size, subset)
%   	dir{}p	- move along {xyz} axis toward point	(step size, subset, point)
%   	dir{}	- move along {xyz} axis			(step size, subset)
%
%   	diff{ef}- diffuse on edges/faces for spacing	(diffusion weight)
%   	edge	- move toward mean of edge neighbors	(step size, % crit dist)
%   	face	- move toward mean of face neighbors	(step size, % crit dist)
%
%   	prune	- remove surface that isn't on target	(pass step)
%   	subf	- subdivide surface (w/ surface freq)
%    	    	  (pass step, % cutoff freq, energy thrshld, resamp factor, max rate,
%   	    	   opt args, xform args)
%   	subp 	- subdivide surface (periodically)	(pass step, samp factor)
%
%   	errf	- compute spectral error (compare to fixed spectrum)
%    	    	  (pass step, fintr, fextr, minerr, xform args)
%
%   	light	- light source position     	    	(position)
%   	surf_e, surf_f, surf_em, surf_fm    	    	(pass step)
%	    	- display or make movie of surface
%   	freq, freq_d, freq_m, freq_dm 	    
%   	    	  (pass step, xform args, comp args, [fintr])
%	    	- display or make movie of spectrum or difference spectrum
%   	    	  (generate fintr if not given)

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - CPU time spent in sf_gen_map, sf_find_neighbors, sf_sub_surf

% - FLOPS is unreliable, since it is not affected by C code

% - ideally surface spectrum would be computed at most once per pass,
%   but different sample sets are used for different purposes...

% - printed output reflects number of points at beginning of pass;
%   this is OK since resampling happens at the end - nothing much has changed

% - step size is controlled by surface spacing (for stability)
%   - explicitly for DIR*
%   - implicitly for DIFF*
%   - ignored for DATA* (instead, fraction of distance to nearby point)

% - DATA* search dist depends on target spacing
% - DIRN wanders dangerously at borders of open surfaces
% - EDGE/FACE transform critical distance (relative to ideal edge length)
%     CD == 1 maintains current spacing (most points don't move)
%     CD <> 1 should slowly enlarge/reduce total area
%   	    	too slow and localized to adjust overall sampling rate
%     CD == 0.9 provides tension for smoothing
%   	    	CD < 1 provides tension for smoothing
% - DIFFE(0.125) is similar to EDGE(0.5,0.9) 


%%% THINGS TO DO

% ? avoid steps without difference spectra by assuming unknown values are zero
%   (or by recomputing spectra after resampling)
%   - current approach may be missing significant energy
%   - rerun simulations?

% ? step size methods with limits based on sampling rate, etc.
% ? time-varying weights (? pairs of weights in dterms)
% ? use FMINS instead of current iteration

% ? DATA* - limit movement to smaller of point distance or step size
%   - may not be necessary since neighboring points will often move together
% - EDGE/FACE
%   ? better transforms (log, x-1, 1-1/x, other combinations)
%   ? if push is proportional to distance, nearby points won't move apart much


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 7) 
    help sf_def_surf; return;
elseif size(intr,2) ~= 2
    error(['INTR must have 2 columns, not ', num2str(size(intr, 2)), '.']);
elseif size(extr, 2) ~= 3
    error(['EXTR must have 3 columns, not ', num2str(size(extr, 2)), '.']);
elseif size(intr, 1) ~= size(extr, 1)
    error(['INTR and EXTR must have equal number of rows (', ...
           num2str(size(intr, 1)), ' ~= ', num2str(size(extr, 1)), ').']);
elseif size(brdr, 2) ~= 1
    error(['BRDR must have 1 columns, not ', num2str(size(brdr, 2)), '.']);
elseif size(intr, 1) ~= size(brdr, 1)
    error(['BRDR and INTR must have equal number of rows (', ...
           num2str(size(brdr, 1)), ' ~= ', num2str(size(intr, 1)), ').']);
elseif size(edge, 2) ~= 2
    error(['EDGE must have 2 columns, not ', num2str(size(edge, 2)), '.']);
elseif size(face, 2) ~= 3
    error(['FACE must have 3 columns, not ', num2str(size(face, 2)), '.']);
elseif all(ndims(target) ~= [2 3])
    error(['TARGET must have 2 or 3 dimensions, not ', num2str(ndims(target)), '.']);
elseif ndims(target) == 2 & size(target, 2) ~= 3
    error(['TARGET must have 2 columns, not ', num2str(size(target, 2)), '.']);
elseif ~iscell(dterms)
    error(['DTERMS must be a cell array.']);
end

% commented this and replaced with pflag=0 to fix error - CLK 2004-02-23
% if strcmp(get(0,'Profile'),'on'), pflag=1; else pflag=0; end
pflag=0;

%%% prepare for deformation
bflag = 0;
oclock = clock;
err = [];
nrby = logical(zeros(size(extr,1), 1));
pass = 0; pass_max = 0;
  pt = 0;   pt_max = 0;
  pp = 0;   pp_max = 0;


subf_err = [0 0 0 0 0 0];
errf_err = 0; 
errf_minerr = 0;
errt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deformation starts here

while 1
    sflag = 0;
    pass = pass + 1;
      pt = size(intr,1);
      pp = pp + pt;

    % reset edge/face info to be recomputed at most once each pass
    % NOTE: these are based on extrinsic coordinates
    evec = []; elen = []; 
    fvec = []; flen = []; farea = [];

    for term = dterms, term = term{1};
    	% echo terms during first pass
%    	if pass==1, disp(term), end
    	
    	switch lower(term{1})

    	    case 'batch', bflag = 1;
    	    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% number of passes, points, passpoints (control term)
    	    case 'pass',
    		if pass==1
    		    sf_def_surf_check_args(term,2); pass_max = term{2}; 
    		end
    	    	if pass >= pass_max, sflag = 1; end
    	    	
    	    case 'point',
    	    	if pass==1
    	    	    sf_def_surf_check_args(term,2); pt_max = term{2};
    	    	end
    	    	if pt >= pt_max, sflag = 1; end
    	    	
    	    case 'passpoint',
    		if pass==1
    		    sf_def_surf_check_args(term,2); pp_max = term{2}; 
    		end
    	    	if pp >= pp_max, sflag = 1; end


    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% step size (control term) (size(ssize) == size(extr))
    	    case 'ssize',
    	    	if pass==1, sf_def_surf_check_args(term,2); end
    	    	if isempty(evec)
    		    [evec,elen] = sf_def_surf_edge_info(extr,edge);
		end

    	    	switch lower(term{2})
    	    	    case 'fixed',   	    	    	%%% fixed step size
    	    	    	if pass==1, sf_def_surf_check_args(term,3); end
    	    	    	ssize = term{3};
    	    	    case 'gmean',   	    	    	%%% global mean edge length
    	    	    	ssize = mean(elen);
    	    	    case 'gmin',    	    	    	%%% global min edge length
    	    	    	ssize = min(elen);
    	    	    case 'lmean',   	    	    	%%% local mean edge length
    		    	ssize = sf_collect(edge(:), elen(:,[1 1]), ...
    		    	    	    	   'full', 'mean', [size(extr,1) 1]);
    	    	    case 'lmin',    	    	    	%%% local min edge length
    		    	ssize = sf_collect(edge(:), elen(:,[1 1]), ...
    		    	    	    	   'full', 'min',  [size(extr,1) 1]);
    	    	end % switch
    	    	
    	    	%%% apply upper limit to step size, if specified
    	    	if length(term) > 2
    	    	    ssize(ssize > term{3}) = term{3};
    	    	end
    	    	%%% expand to appropriately sized matrix
    	    	if prod(size(ssize)) == 1,
    	    	    ssize = ssize(ones(size(extr)));
    	    	elseif all(size(ssize) == [size(extr,1) 1])
    	    	    ssize = ssize(:,[1 1 1]);
    	    	end
    	    	    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% scattered data method
    	    %%% - set nearby surface points, move toward nearest target point
    	    case 'datas',
    	    	if pass==1,
    	    	    sf_def_surf_check_args(term,4);
    		    if ndims(target) ~= 2
    	    		error(['Term ', term{1}, ' requires scattered target data.']);
    		    end
    		end

    	    	% find number and mean position of target points near surface
    	    	% flag nearby points, move surface toward mean of nearby targets
    	    	% NOTE: this is very time-consuming
    	    	[datas_i,datas_d,datas_c,datas_m] ...
    	    	    = sf_find_neighbors(extr, target, term{4}, term{3});
                % only want to affect points that have neighbors
                datas_c = datas_c > 0;
                nrby(datas_c) = 1;
    	    	% step size depends on distance, so don't normalize
    	    	extr(datas_c,:) = extr(datas_c,:) + ...
    	    	    	    	  term{2} * (datas_m(datas_c,:) - extr(datas_c,:));

    	    	% compute target error (for validation)
    	    	errt = mean(datas_d);


    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% volume data method
    	    %%% - set nearby surface points, move toward mean of nearest target points
    	    case 'datav'
    	    	%%% check and/or recompute as needed
    	    	if pass==1,
    		    sf_def_surf_check_args(term,3);
    		    if ndims(target) ~= 3
    	    		error(['Term ', term{1}, ' requires volumetric target data.']);
    		    end
    		    tmin = []; tmax = [];
    	    	end

    	    	if size(tmin,1) ~= size(extr,1)
   		    tmin =      [1 1 1] + term{3}; tmin = tmin(ones(size(extr,1),1),:);
   		    tmax = size(target) - term{3}; tmax = tmax(ones(size(extr,1),1),:);
    	    	end

    	    	% EINS = indices of surface points inside target volume
    	    	eins = find(all((extr > tmin) & (extr < tmax), 2));
    	    	
    	    	%%% if inside points, search range of neighbors for target
    	    	if size(eins,1)
    	    	    % EINT = 3D coordinates, EIND* = 1D coordinates
    	    	    eint = zeros(size(eins,1), 3);
    	    	    % NCNT = neighbor counts, NPOS = coordinate sums (for inside points)
    	    	    ncnt = zeros(size(eins,1), 1);
    	    	    npos = zeros(size(eins,1), 3);
    	    	    extr1 = extr(eins,1); extr2 = extr(eins,2); extr3 = extr(eins,3);
    	    	    range = -term{3}:term{3};

    	    	    for j3 = range
    	    	    	eint(:,3) = round(extr3 + j3);
    	    	    	eind3 = (eint(:,3) - 1) * size(target,2);
    	    	    	for j2 = range
    	    	    	    eint(:,2) = round(extr2 + j2);
    	    	    	    eind2 = (eint(:,2) - 1 + eind3) * size(target,1);
    	    	    	    for j1 = range
    	    	    	    	eint(:,1) = round(extr1 + j1);
    	    	    	    	eind1 = (eint(:,1) + eind2);
    	    	    	    	% EIND = flag for target neighbors of surface points
    	    	    	    	eind = target(eind1);	    	    	
    	    	    		% add coordinates to npos, increment ncnt
    	    	    		npos(eind,:) = npos(eind,:) + eint(eind,:);
    	    	    		ncnt(eind)   = ncnt(eind)   + 1;
    	    	    	    end % for J3
    	    		end % for J2
    	    	    end % for J1
                    
    	    	    % EINS = indices of surface points with neighbors
                    ncntpos = ncnt>0;
    	    	    eins = eins(ncntpos);
    	    	    % flag nearby surface points (those with neighbors)
    	    	    nrby(eins) = 1;
    	    	    % NPOS = mean position of neighbors for each surface point
    	    	    npos = npos(ncntpos,:) ./ ncnt(ncntpos,[1 1 1]);
    	    	    % don't normalize, since step size depends on distance
    	    	    extr(eins,:) = extr(eins,:) + term{2} * (npos - extr(eins,:));

    	    	    % compute target error (for validation)
    	    	    errt = (sum(sqrt(sum((npos - extr(eins,:)).^2, 2))) + ...
    	    	    	    sum(2*term{3}*(size(extr,1) - size(eins,1)))) / size(extr,1);
    	    	else
    	    	    errt = 2*term{3};
    	    	end % if inside points
   	    	    

    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	    %%% direction methods - move subset of points in specified direction

    	    case 'dirr', 
    	    	if pass==1, 
    	    	    sf_def_surf_check_args(term,3); 
    	    	    if exist('ssize') ~= 1, error(['Term ', term{1}, ' requires step size.']); end
    	    	end
    	    	dir_pnts = sf_def_surf_dir_pnts(term{3},nrby);
   	    	extr(dir_pnts,:) = extr(dir_pnts,:) + ...
    	    	    term{2} * ssize(dir_pnts,:) .* randn(size(extr(dir_pnts,:)));

    	    case {'dirx','diry','dirz'},
    	    	if pass==1, 
    	    	    sf_def_surf_check_args(term,3); 
    	    	    if exist('ssize') ~= 1, error(['Term ', term{1}, ' requires step size.']); end
    	    	end
    	    	switch lower(term{1})
    	    	    case 'dirx', dir_c = 1;
    	    	    case 'diry', dir_c = 2;
    	    	    case 'dirz', dir_c = 3;
    	    	end
    	    	dir_pnts = sf_def_surf_dir_pnts(term{3},nrby);
    	    	extr(dir_pnts,dir_c) = extr(dir_pnts,dir_c) + ...
    	    	    term{2} * ssize(dir_pnts,dir_c);

    	    case 'dirn',
    	    	if pass==1, 
    	    	    sf_def_surf_check_args(term,3);
    	    	    if exist('ssize') ~= 1, error(['Term ', term{1}, ' requires step size.']); end
    	    	end
    	    	dir_pnts = sf_def_surf_dir_pnts(term{3},nrby);
    	    	fnrm = cross(extr(face(:,1),:) - extr(face(:,2),:), ...
    	    	    	     extr(face(:,2),:) - extr(face(:,3),:));
    	    	dirn = [ sf_collect(face(:), fnrm(:,[1 1 1]), 'full', 'sum', [size(extr,1),1]), ...
    	    	    	 sf_collect(face(:), fnrm(:,[2 2 2]), 'full', 'sum', [size(extr,1),1]), ...
    	    	    	 sf_collect(face(:), fnrm(:,[3 3 3]), 'full', 'sum', [size(extr,1),1]) ];
    	    	extr(dir_pnts,:) = extr(dir_pnts,:) + ...
    	    	    term{2} * ssize(dir_pnts,:) .* ...
    	    	      sf_def_surf_norm(dirn(dir_pnts,:));
    	    	
    	    case {'dirxp','diryp','dirzp','dirxyp','dirxzp','diryzp','dirxyzp','dirp'},
    	    	if pass==1,
    	    	    sf_def_surf_check_args(term,4);
    	    	    if exist('ssize') ~= 1, error(['Term ', term{1}, ' requires step size.']); end
    	    	    dir_p = [];
    	    	end
    	    	switch lower(term{1})
    	    	    case 'dirxp',   	    	dir_c = [1    ];
    	    	    case 'diryp',   	    	dir_c = [  2  ];
    	    	    case 'dirzp',   	    	dir_c = [    3];
    	    	    case 'dirxyp',  	    	dir_c = [1 2  ];
    	    	    case 'dirxzp',  	    	dir_c = [1   3];
    	    	    case 'diryzp',  	    	dir_c = [  2 3];
    	    	    case {'dirxyzp','dirp'},	dir_c = [1 2 3];
    	    	end
    	    	dir_pnts = sf_def_surf_dir_pnts(term{3},nrby);
    	    	if size(dir_p,1) ~= size(extr,1)
    	    	    dir_p = term{4};
    	    	    dir_p = dir_p(ones(size(extr,1),1),:);
    	    	end
    	    	extr(dir_pnts,dir_c) = extr(dir_pnts,dir_c) + ...
    	    	    term{2} * ssize(dir_pnts,dir_c) .* ...
    	    	      sf_def_surf_norm(dir_p(dir_pnts,dir_c) - extr(dir_pnts,dir_c));

    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	    %%% diffusion methods - balance spacing between surface points
   	    %%%   distant points (longer edges, larger areas) have more weight
   	    case 'diffe',
    	    	%%% check and/or recompute as needed
    	    	if pass==1 sf_def_surf_check_args(term,2); diffe_ndx = []; end
    	    	if isempty(evec)
    	    	    [evec,elen] = sf_def_surf_edge_info(extr,edge);
    	    	end
    	    	if size(diffe_ndx,1) ~= size(extr,1)
    		    diffe_ndx = 1:size(extr,1);
    		    diffe_row = reshape(edge(:,[1 2]), prod(size(edge)), 1);
    		    diffe_col = reshape(edge(:,[2 1]), prod(size(edge)), 1);
    	    	end
    	    	
    	    	%%% EDIFF = edge-based diffusion matrix (opposite weights from SF_OPT_INTR)
    	    	%%% - each edge length counts toward both end points
    	    	%%% - border points only on diagonals
    	    	ediff = sparse(diffe_ndx,diffe_ndx,  brdr) + ...
    	    	        sparse(diffe_ndx,diffe_ndx, ~brdr) * ...
    	    	    	    sparse(diffe_row,diffe_col, elen(:,[1 1]), ...
    	    	    	   	   size(extr,1),size(extr,1));
    	    	% normalize (row sums == 1) and weight diffusion matrix
    	    	ediff = sparse(diffe_ndx,diffe_ndx, 1 ./ sum(ediff,2)) * ediff;
    	    	ediff = (1 - term{2}) * sparse(diffe_ndx,diffe_ndx, 1) + ...
    	    	    	(    term{2}) * ediff;
    	    	extr = ediff * extr;    	    	    	
   	    

   	    case 'difff',
    	    	%%% check and/or recompute as needed
    	    	if pass==1 sf_def_surf_check_args(term,2); difff_ndx = []; end
    	    	if isempty(fvec)
    	    	    [fvec,flen,farea] = sf_def_surf_face_info(extr,face);
    	    	end
    	    	if size(difff_ndx,1) ~= size(extr,1)
    		    difff_ndx = 1:size(extr,1);
    		    difff_row = reshape(face(:, [1 1 2 2 3 3]), 2*prod(size(face)), 1);
    		    difff_col = reshape(face(:, [2 3 1 3 1 2]), 2*prod(size(face)), 1);
    	    	end
    	    	
       	    	%%% FDIFF = face-based diffusion matrix (opposite weights from SF_OPT_INTR)
       	    	%%% - each face area counts toward all three corner points
       	    	%%% - border points only on diagonal
    	    	fdiff = sparse(difff_ndx,difff_ndx,  brdr) + ...
    	    	        sparse(difff_ndx,difff_ndx, ~brdr) * ...
    	    	    	    sparse(difff_row,difff_col, farea(:,[1 1 1 1 1 1]), ...
    	    	    	    	   size(extr,1),size(extr,1));
    	    	% normalize (row sums == 1) and weight diffusion matrix
    	    	fdiff = sparse(difff_ndx,difff_ndx, 1 ./ sum(fdiff,2)) * fdiff;
    	    	fdiff = (1 - term{2}) * sparse(difff_ndx,difff_ndx, 1) + ...
    	    	        (    term{2}) * fdiff;
    	    	extr = fdiff * extr;
   	    

    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% edge method - attract/repel edge neighbors
    	    case 'edge',
    	    	%%% check and/or recompute as needed
    	    	if pass==1, sf_def_surf_check_args(term,3); end
    	    	if isempty(evec)
    	    	    [evec,elen] = sf_def_surf_edge_info(extr,edge);
    	    	end

    	    	% EWGT = edge weights (-1...0 to repel, 0...1 to attract)
    	    	ewgt = elen / (term{3} * median(elen));
    	    	ewgt(ewgt < 1) = ewgt(ewgt < 1) - 1;    	    	
    	    	ewgt(ewgt > 1) = 1 - (1 ./ ewgt(ewgt > 1));
    	    	% EDIR = vector to mean of edge neighbors
    	    	edir = evec .* ewgt(:,[1 1 1]);
    	    	edir = [edir; -edir];
    	    	edir = [sf_collect(edge(:), edir(:,1), 'full', 'mean', [size(extr,1),1]), ...
    	    	    	sf_collect(edge(:), edir(:,2), 'full', 'mean', [size(extr,1),1]), ...
    	    	    	sf_collect(edge(:), edir(:,3), 'full', 'mean', [size(extr,1),1]) ];
    	    	% step size depends on weighting, so don't normalize
    	    	extr(~brdr,:) = extr(~brdr,:) + term{2} * edir(~brdr,:);


    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% face method - attract/repel face neighbors
    	    case 'face',
    	    	%%% check and/or recompute as needed
    	    	if pass==1, sf_def_surf_check_args(term,3); end
    	    	if isempty(fvec)
    	    	    [fvec,flen,farea] = sf_def_surf_face_info(extr,face);
    	    	end
    	    	
    	    	% FWGT = face weights (-1...0 to repel, 0...1 to attract)
    	    	fwgt = farea / (term{3} * median(farea));
    	    	fwgt(fwgt < 1) = fwgt(fwgt < 1) - 1;
    	    	fwgt(fwgt > 1) = 1 - (1 ./ fwgt(fwgt > 1));
    	    	% FDIR = vector to mean of face neighbors
    	    	fdir  = reshape(fvec .* fwgt(:,[1 1 1],[1 1 1]), 3*size(face,1), 3);
    	    	fdir = [ sf_collect(face(:), fdir(:,1), 'full', 'mean', [size(extr,1),1]), ...
    	    	    	 sf_collect(face(:), fdir(:,2), 'full', 'mean', [size(extr,1),1]), ...
    	    	    	 sf_collect(face(:), fdir(:,3), 'full', 'mean', [size(extr,1),1]) ];
    	    	% step size depends on weighting, so don't normalize
    	    	extr(~brdr,:) = extr(~brdr,:) + term{2} * fdir(~brdr,:);
    	    	

    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% prune method
    	    case 'prune',
    	    	if pass==1, sf_def_surf_check_args(term,2); end
    	    	if rem(pass,term{2}) == 0
    	    	    % compute cumulative sum to get new indices
    	    	    csnn = cumsum(nrby);
    	    	    edge(sf_find_vals(find(~nrby),edge),:) = [];
    	    	    face(sf_find_vals(find(~nrby),face),:) = [];
    	    	    edge = csnn(edge);
    	    	    face = csnn(face);
    	    	    intr = intr(nrby,:); 
    	    	    extr = extr(nrby,:); 
    	    	    brdr = brdr(nrby,:);
   	    	    nrby = nrby(nrby,:);
    	    	end

    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% frequency-based subdivision method
    	    case 'subf', 
    	    	if pass==1
    	    	    sf_def_surf_check_args(term,8); 
    	    	    subf_cutoff = term{3};
    	    	    subf_thresh = term{4};
    	    	    subf_resamp = term{5};
    	    	    subf_maxsam = term{6};
    	    	      subf_opt_args = term{7};
    	    	    subf_xform_args = term{8};
    	    	    subf_area = []; subf_fintr = []; subf_farea = []; 
    	    	    subf_fextr0 = [];
    	    	end
    	    	
    	    	if rem(pass,term{2}) == 0
    	    	    %%% recompute freq samples after resampling
    	    	    if size(subf_fintr,1) < size(intr,1)
    	    	    	subf_fintr = sf_gen_intr(ceil(sqrt(size(intr,1))), 'rectq', ...
    	    	    	    	    	    	 ceil(sqrt(size(intr,1))), 0.01);
    	    	    	subf_farea = sparse(1:size(subf_fintr,1),1:size(subf_fintr,1), ...
    	    	    	    	    	    sf_gen_area(subf_fintr, size(subf_fintr,1), 'con'));
    	    	    end
    	    	    
    	    	    %%% optimize coordinates and get surface frequency
    	    	    intr = sf_opt_intr(subf_opt_args{:}, intr,extr,brdr,edge,face);
    	    	    subf_area = sparse(1:size(intr,1),1:size(intr,1), ...
    	    	    	    	       sf_gen_area(intr,1,'face',face));
    	    	    [map, fmap] = sf_gen_map(intr,subf_fintr,'ft');
    	    	     map =  map * subf_area;
    	    	    fmap = fmap * subf_farea;
    	    	    subf_fextr = sf_xform(sf_gxform(extr,'nrm',[]), ...
    	    	    	    		 map, fmap, subf_xform_args{:});

    	    	    %%% reset comparison spectrum and cumulative sum after resampling
 	    	    if size(subf_fextr0,1) ~= size(subf_fextr,1),
 	    	    	subf_fextr0 = subf_fextr;
    	    	    	subf_csum   = 0;
 	    	    end 

    	    	    %%% get composite spectrum of successive difference
    	    	    [subf_cintr, subf_cextr] = sf_comp(subf_fintr, subf_fextr - subf_fextr0, ...
    	    	    	    	    	    	       '1dm','pow','bin');
    	    	    subf_fextr0 = subf_fextr;
    	    	    
    	    	    %%% add high frequency energy to cumulative sum  
    	    	    fc = subf_cutoff * subf_cintr(end);
    	    	    eda = sum(subf_cextr) + eps;
    	    	    edh = sum(subf_cextr(subf_cintr >= fc));
    	    	    subf_csum = subf_csum + edh;
    	    	    subf_err = [edh, eda, edh/eda, fc, subf_csum, subf_thresh];
    	    	    	    
    	    	    %%% if high freq energy exceeds threshold, subdivide surface
    	    	    if (subf_maxsam > length(intr)) & (subf_csum >= subf_thresh)
%    	    		[intr,extr,brdr,edge,face] = ...
%    	    	    	    sf_sub_surf(intr,extr,brdr,edge,face,'pfac',subf_resamp);
    	    		[intr,extr,brdr,edge,face] = ...
    	    	    	    sf_sub_surf(intr,extr,brdr,edge,face,'ptot', ...
    	    	    	    	    	min(subf_maxsam,subf_resamp*length(intr)));
    	    		nrby(end:size(intr,1)) = 0;
    	    	    end % if energy > threshold    	    	
    	    	end


    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% periodic subdivision method
    	    case 'subp', 
    	    	if pass==1, sf_def_surf_check_args(term,3); end
    	    	if rem(pass,term{2}) == 0
    	    	    [intr,extr,brdr,edge,face] = ...
    	    	    	sf_sub_surf(intr,extr,brdr,edge,face,'pfac',term{3});
    	    	    nrby(end:size(intr,1)) = 0;
    	    	end
    	    	
    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% spectral error method
    	    case 'errf',
    	    	if pass==1
    		    sf_def_surf_check_args(term,6);
    		    errf_fintr      = term{3};
    		    errf_fextr0     = term{4};
    		    errf_minerr     = term{5};
    		    errf_xform_args = term{6};
    		    errf_farea = sparse(1:size(errf_fintr,1),1:size(errf_fintr,1), ...
    		    	    	    	sf_gen_area(errf_fintr, size(errf_fintr,1), 'con'));
    	    	end
    	    
    	    	if pass==1 | rem(pass,term{2})==0
    	    	    errf_area = sparse(1:size(intr,1),1:size(intr,1), ...
    	    	    	    	       sf_gen_area(intr,1, 'face',face));
    	    	    [map, fmap] = sf_gen_map(intr,errf_fintr, 'ft');
    	    	     map =  map * errf_area;
    	    	    fmap = fmap * errf_farea;
    	    	    errf_fextr = sf_xform(sf_gxform(extr,'nrm',[]), ...
    	    	    	    	    	   map,fmap, errf_xform_args{:});
    	    	    errf_err = sum(sf_norm(errf_fextr0 - errf_fextr,'energy',errf_farea));
    	    	    if errf_err < errf_minerr, sflag = 1; end
    	    	    
    	    	end
    	    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% surface visualization methods
       	    case {'surf_e','surf_f','surf_em','surf_fm'},
    		if pflag, profile off; end
    	    	if pass==1,
    	    	    %%% initialize variables, etc.
    		    sf_def_surf_check_args(term,2); 
    		    surf_fh1 = figure('Name',			'sf_def_surf', ...
    	    			      'PaperPositionMode',	'auto', ...
    	    			      'Position',		[0 400 500 500], ...
    	    			      'DefaultLineLineWidth',	0.25, ...
    	    			      'DefaultLineMarkerSize',	2);
    	    	    surf_ph1 = [];
    		    view(-45,30); rotate3d; sf_tool all

    		    %%% plot bounding box and 500 target points for reference
    		    if ndims(target) == 2
    		    	bbox = num2cell([ min(target) max(target) ]);
    	    		plotpts = sf_sample(target,500);
    	    		ppx = plotpts(:,1); ppy = plotpts(:,2); ppz = plotpts(:,3);
    		    elseif ndims(target) == 3
    		    	bbox = num2cell([ 1 1 1 size(target) ]);
    	    		plotpts = sf_sample(find(target), 500);
    	    		[ppx,ppy,ppz] = ind2sub(size(target),plotpts);
    		    end
    	    	    hold on
    	    	    sf_def_surf_plot_box(bbox{:});
    		    plot3(ppx, ppy, ppz, 'g+');
    		    hold off
    		    surf_ah1 = gca;

    		    aran = [ min(axis) max(axis) ];
    		    axis([ aran aran aran ]); axis('equal', 'vis3d');
    		    zoom fill;
    		    xlabel('x'); ylabel('y'); zlabel('z');
    		    switch term{1}
    	    		case {'surf_e','surf_f'},
    	    		case {'surf_em','surf_fm'},
    	    		    set(surf_ah1, 'Units', 'pixels');
    	    		    smovie_rect = get(surf_ah1, 'Position');
    	    		    smovie      = moviein(floor(pass_max/term{2}) + 1, ...
    	    		    	    	    	   surf_fh1, smovie_rect);
    	    		    smovie(:,1) = getframe(surf_fh1, smovie_rect);
    	    	    end % switch
    	    	end % if pass==1

    	    	if pass==1 | rem(pass,term{2})==0
   		    switch term{1}
    	    		case {'surf_e','surf_em'}, ecolor =  'flat'; fcolor = 'none';
    	    		case {'surf_f','surf_fm'}, ecolor = [0 0 0]; fcolor = 'flat';
    	    	    end
    	    	    delete(surf_ph1);
     		    surf_ph1 = patch('Faces',		face, ...
    				     'Vertices',	extr, ...
    				     'FaceVertexCData',	double(nrby), ...
    				     'EdgeColor',	ecolor, ...
    				     'FaceColor',     	fcolor, ...
    				     'BackFaceLighting','reverselit', ...
    				     'EdgeLighting',	'none', ...
    				     'FaceLighting',	'none', ...
    				     'Parent',	    	surf_ah1, ...
    				     'Clipping',      	'off');

    		    switch term{1}
    	    		case {'surf_e','surf_f'},
    	    	    	    pause(0.01);
    			case {'surf_em','surf_fm'},
    	    	    	    smovie(:,floor(pass/term{2}) + 1) = ...
    	    	    	    	getframe(surf_fh1, smovie_rect);
    	    	    end
    	    	end % if 
    	    	if pflag, profile on; end

    	    	
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% light method
    	    case 'light',
    	    	if pass==1
    		    sf_def_surf_check_args(term,2);
    		    if exist('surf_fh1') == 1, figure(surf_fh1); end
    	    	    light('Position',term{2});
    	    	end
    	    
    	    
    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    %%% surface frequency visualization methods
    	    case {'freq', 'freq_m', 'freq_d', 'freq_dm'},
    		if pflag, profile off; end
    	    	if pass==1,
    	    	    %%% initialize variables, etc.
    		    freq_fh1 = [];
    		    freq_area = []; freq_fintr = []; freq_farea = [];
    		    freq_fextr0 = []; freq_fextr = [];
    		    sf_def_surf_check_args(term,4);
    		    freq_xform_args = term{3};
    		     freq_comp_args = term{4};
    		    if length(term) >= 5
    		    	freq_fintr = term{5};
    			if size(freq_fintr,2) ~= 2
    		    	    error(['Term ', term{1}, ' has invalid argument.']);
    			end
    		    end
    		end % if pass==1

    	    	if pass==1 | rem(pass,term{2})==0
    	    	    %%% recompute freq samples if not given and intr has changed
    	    	    if length(term) < 5 & (size(freq_fintr,1) < size(intr,1))
    		    	freq_fintr = sf_gen_intr(ceil(sqrt(size(intr,1))), 'rectq', ...
    		    	    	    	    	 ceil(sqrt(size(intr,1))), 0.01);
    		    end
    		    %%% recompute freq areas if freq samples have changes
    		    if size(freq_farea,1) ~= size(freq_fintr,1)
    			freq_farea = sparse(1:size(freq_fintr,1),1:size(freq_fintr,1), ...
    		    	    		    sf_gen_area(freq_fintr, size(freq_fintr,1), 'con'));
    	    	    end

    	    	    %%% compute and plot surface frequency
    		    freq_area = sparse(1:size(intr,1),1:size(intr,1), ...
    		    	    	       sf_gen_area(intr,1,'face',face));
    	    	    [map, fmap] = sf_gen_map(intr,freq_fintr,'ft');
    	    	     map =  map * freq_area;
    	    	    fmap = fmap * freq_farea;
    	    	    freq_fextr = sf_xform(sf_gxform(extr,'nrm',[]), ...
    	    	    	    	          map,fmap, freq_xform_args{:});
    	    	    switch term{1}
    	    	    	case {'freq'  ,'freq_m' },
    	    	    	case {'freq_d','freq_dm'}, 
    	    	    	    if size(freq_fextr0,1) ~= size(freq_fextr,1),
    	    	    	    	freq_fextr0 = freq_fextr;
    	    	    	    end 
    	    	    	    freq_fextr = freq_fextr - freq_fextr0;
    	    	    end

    	    	    % sf_comp is mostly graphics calls, so ignore while profiling
    	    	    freq_fh1 = sf_comp(freq_fintr, freq_fextr, freq_fh1, freq_comp_args{:});

    		    switch term{1}
    	    		case {'freq', 'freq_d'},
    	    	    	    pause(0.01);
    			case {'freq_m', 'freq_dm'},
    			    if pass==1
    	    			set(freq_fh1, 'Units', 'pixels');
    	    			fmovie_rect = get(freq_fh1, 'Position');
    	    			fmovie_rect(1:2) = [0 0];
    	    			fmovie      = moviein(floor(pass_max/term{2}) + 1, ...
    	    		    	    	    	      freq_fh1, fmovie_rect);
    	    		    	fmovie(:,1) = getframe(freq_fh1, fmovie_rect);
    			    end 
    	    	    	    fmovie(:,floor(pass/term{2}) + 1) = ...
    	    	    	    	getframe(freq_fh1, fmovie_rect);
    	    	    end % switch
    	    	end % if
    	    	if pflag, profile on; end

    	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	    otherwise,
    		if pass==1, disp(['Unknown term: ', term{1}, '.']); end

    	    end % switch TERM

    end % for TERM = DTERMS


    %%% save error values, print status messages
    format0 = '\n  secs    pass       nrby pt          pp        errt        errf         high    tot     frac    cutoff     cum    thresh \n';
    format1 = '%7.2f  %3d/%3d  %4d/%4d/%4d  %5d/%5d  %6.4f  %6.4f/%6.4f  %8.6f/%8.6f (%5.3f) > %5.2f  %8.6f/%8.6f \n';
    
    format = '%6.2f sec, %3d/%3d pass, %4d/%4d/%4d nrby pt, %5d/%5d pp, %f errt, %f/%f errf, %f/%f (%5.3f) df > %5.2f, %f/%f cum \n';

    errnew = [etime(clock,oclock), ...
    	      pass, pass_max, size(find(nrby),1), pt, pt_max, pp, pp_max, ...
    	      errt, errf_err, errf_minerr, subf_err ];
    err    = [err; errnew];
    
    if 1 | ~bflag
    	if pass==1, fprintf(1, format0); end
    	fprintf(1,format1, err(pass,:));
    end  

    %%% if sflag was set, play movies, adjust limits, etc.
    if sflag == 1
    	if bflag, return; end
    	if pflag, profile off; end
    	clk_beep(2);
    	while sflag == 1
    	    if exist('fmovie'),   fprintf(1,'<F>-movie [n fps], '); end
    	    if exist('smovie'),   fprintf(1,'<S>-movie [n fps], '); end
    	    if pass_max > 0, fprintf(1,'<P> [passes], ');      end
    	    if	 pt_max > 0, fprintf(1,'<N> [points], ');      end
    	    if   pp_max > 0, fprintf(1,'<PN> [passpoints], '); end
    	    fprintf(1,'<C>ont, <D>ebug, <Q>uit ');
    	    [ans1,ans2] = strtok(input('? ','s')); 
    	    if isempty(ans1), ans1 = ''; end
    	    if isstr(ans2), ans2 = str2num(ans2); end
    	    switch upper(ans1)
    		case 'F',
    		    if exist('fmovie')
    		    	if ~(isnumeric(ans2) & prod(size(ans2)) == 1) ans2 = 2; end
    			movie(freq_fh1, fmovie, ans2, 5, fmovie_rect);
    		    end
    		case 'S', 
    		    if exist('smovie')
    		    	if ~(isnumeric(ans2) & prod(size(ans2)) == 1) ans2 = 2; end
    			movie(surf_fh1, smovie, ans2, 5, smovie_rect);
    		    end

    		case 'P',
    	    	    if isnumeric(ans2) & prod(size(ans2)) == 1
    	    		pass_max = pass_max + ans2;
    	    	    end
    	    	    fprintf(1,'%3d/%3d passes \n', pass, pass_max);
    	    	case 'N',
    	    	    if isnumeric(ans2) & prod(size(ans2)) == 1
    	    	    	pt_max = pt_max + ans2;
    	    	    end
    	    	    fprintf(1,'%4d/%4d points \n', pt, pt_max);
    		case {'NP','PN','PP'},
    	    	    if isnumeric(ans2) & prod(size(ans2)) == 1
    	    		pp_max = pp_max + ans2;
    	    	    end
    	    	    fprintf(1,'%3d/%3d passpoints \n', pp, pp_max);

    		case  'C',   	break;
    		case  'D',   	keyboard;
    		case {'Q','R'}, return;
    	    end % switch ans1
    	end % while sflag
    	if pflag, profile on; end
    end % if sflag

end % while 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions start here

function sf_def_surf_check_args(term,num)
    if length(term) < num
    	error(['Term ', term{1}, ' requires ', num2str(num), ' argument(s).']);
    end

function dir_pnts = sf_def_surf_dir_pnts(term,nrby);
    switch lower(term)
    	case 'all',         	dir_pnts =  ones(size(nrby));
    	case {'dstnt','nnrby'}, dir_pnts =           ~nrby;
    	case {'ndstnt','nrby'}, dir_pnts =            nrby;
    	case 'none',        	dir_pnts = zeros(size(nrby));
    	otherwise,          	error(['Invalid direction mode.']);
    end    	    	    

function [evec,elen] = sf_def_surf_edge_info(pts,edge)
    epts = reshape(pts(edge,:),[size(edge) size(pts,2)]);
    evec = squeeze(diff(epts,1,2));
    % ELEN = length of all edges
    elen = sqrt(sum(evec.^2,2));

function [fvec,flen,farea] = sf_def_surf_face_info(pts,face)
    fpts  = reshape(pts(face,:), [size(face) size(pts,2)]);
    fmean = mean(fpts,2);
    % FVEC = vectors from center to corners of face
    fvec  = fmean(:,[1 1 1],:) - fpts;
    % FMEAN = face centers
    fmean = squeeze(fmean);
    fpts  = reshape(pts(face(:,[1 2 3 1]),:), [size(face,1) 4 size(pts,2)]);
    flen  = sqrt(sum(diff(fpts,1,2).^2,3));
    fs	  = sum(flen,2)/2;
    % FAREA = face areas
       fa = sqrt(prod([fs, fs(:,[1 1 1]) - flen],2));
    farea = sf_collect(face(:), [fa;fa;fa], 'sum', 'full', [size(pts,1) 1]);
 
function v = sf_def_surf_norm(v)
    vn = sqrt(sum(v.^2,2));
    v  = v ./ vn(:,ones(size(v,2),1));
      
function sf_def_surf_plot_box(x1,y1,z1,x2,y2,z2)
    plot3([x1 x2 x2 x1 x1; x1 x2 x2 x1 x1; x1 x1 x1 x1 x1; x2 x2 x2 x2 x2]', ...
    	  [y1 y1 y2 y2 y1; y1 y1 y2 y2 y1; y1 y2 y2 y1 y1; y1 y2 y2 y1 y1]', ...
    	  [z1 z1 z1 z1 z1; z2 z2 z2 z2 z2; z1 z1 z2 z2 z1; z1 z1 z2 z2 z1]', 'g-');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

