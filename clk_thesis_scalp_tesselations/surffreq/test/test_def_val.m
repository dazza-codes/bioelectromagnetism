
% test_def_val.m    - validate SF-based resampling in sf_def_model

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - to get comparable errors, surfaces with different sampling rates
%   need comparable intrinsic ranges.  Since sf_gen_intr adjusts
%   the intr range based on sampling rate, all surfaces must be 
%   generated at the same sampling rate, and subsampled to desired rate

% - to get comparable errors, deformed surfaces need comparable positions.
%   Therefore, pick number of baseline passes for useful model, 
%   and run fixed and resampling models repeatedly past similar point.
%   Collect output statistics from pass with minimum spectral error,
%   (when the model was closest to baseline) and compute averages, etc.
%   Finally, rerun with average number of passes for typical output data

% - using best methods from ANOVA
%   - using DIFC for optimization
%   - using FACE for area weights

%%% RESULTS

% subf 1 0.707 0.000250 2 npoints

%   surf=sheet reps=10, freq: n=20, bas np=40-30, fix np=20-12, res np=10-9 
%   fixed:  10 reps, 12 passes, pp = 4800(   0), err = 0.0170(0.0013) 
%   resam:  10 reps,  9 passes, pp = 2520( 140), err = 0.0208(0.0026) 
%   Prob of significantly different errors (t-test): 1.00

%   surf=sheet reps=02, freq: n=30, bas np=40-30, fix np=30-15, res np=10-10 
%   fixed:   2 reps, 15 passes, pp = 18900(   0), err = 0.0138(0.0013) 
%   resam:   2 reps, 10 passes, pp =  6650( 636), err = 0.0170(0.0004) 
%   Prob of significantly different errors (t-test): 0.84 

% rdir = -0.05
%   surf=sheet reps=05, freq: n=30, bas np=40-30, fix np=30-25, res np=10-20 
%   fixed:   5 reps, 25 passes, pp = 19080( 753), err = 0.0115(0.0012) 
%   resam:   5 reps, 20 passes, pp =  5760( 518), err = 0.0127(0.0022) 
%   Prob of significantly different errors (t-test): 0.66 

% rdir = -0.1
%   resam:   5 reps, 20 passes, pp =  5720( 415), err = 0.0119(0.0009) 
% rdir = -0.15
%   resam:   5 reps, 20 passes, pp =  5460( 365), err = 0.0122(0.0014)


%%% THINGS TO DO
% ? vary resampling factor - probably less meaningful than other parameters
% - rerun simulations with more passes, smaller steps, different range
%   - someday when there is lots of time or very fast hardware

% ? for final run, set bpoints = 60
% - may not want to use sf_nspec for output spectra (WHY OR WHY NOT?)
% - summarize time, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% declarations start here

clear
more off

 plot_flag = 1;  	% boolean flag to plot surfaces
xplot_flag = 1;  	% boolean flag to plot spectra

reps	= 5;

bpass	= 40;	% passes in baseline model
npass	= 30;	% passes in fixed model
rpass   = 30;	% passes in resampling model

bpoints = 40;	% points in baseline model
fpoints = 30;	% points in comparison spectra
npoints = 30;	% points in fixed model
rpoints = 10;	% points in resampling model

cutofflist = 0.400:0.100:0.900;
threshlist = 0.000050:0.000050:0.000250;

  opt_args = {'difcq','deld',0.001};	    % arguments for sf_opt_intr
xform_args = {'traq', 'delf',0.01, 0.1};    % arguments for sf_xform
xform_args = {'dirq'};
 comp_args = {'xyzp','2d'}; 	    	    % arguments for sf_comp

data	= 'mri';
surf	= 'sheet';

format0 = '\n  secs    pass       nrby pt          pp        errt        errf         high    tot     frac    cutoff     cum    thresh \n';
format1 = '%7.2f  %3d/%3d  %4d/%4d/%4d  %5d/%5d  %6.4f  %6.4f/%6.4f  %8.6f/%8.6f (%5.3f) > %5.2f  %8.6f/%8.6f \n';

fid = 1;
fid = fopen('test_def_val.out', 'a');

fprintf(fid, '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid, '%% %4d/%02d/%02d %02d:%02d:%02d ', fix(clock)); 
fprintf(fid, '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');

fprintf(fid,'surf=%s reps=%d, freq: n=%d, base: n=%d p=%d, fixed: n=%d p=%d, resam: n=%d p=%d \n', ...
    	surf, reps, fpoints, bpoints, bpass, npoints, npass, rpoints, rpass);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% selections start here

%%% read data_vals
switch data
    case 'byu',
    	data_vals = sf_read_byu('head.byu');
	data_vals = sf_gxform(data_vals, ...
    	    		     'te',    -(max(data_vals)+min(data_vals))/2, ...
    	    		     'sa',1/max(max(data_vals)-min(data_vals)), ...
    	    		     'rx',pi/2);
    case 'mri',
    	load clkmri003
    	data_vals = sf_find_3d(clkmri003,80,'vol');
    	clear clkmri003
end % switch data


%%% set gxform args and direction terms for various surfaces
switch surf
    case 'cylinder',
    	switch data
    	    case 'byu',
    	    	gxform_args = {'rz',pi/2, 'se',[sqrt(2) sqrt(2), 1]};
    	    case 'mri',
    	end % switch data

    	   dir_term = {'dirxyp', 0.5, 'dstnt', [0 0 0]};
    	  rdir_term = {'dirxyp',-0.1, 'nrby' , [0 0 0]};
    	
    case 'sheet',
    	switch data
    	    case 'byu',
    	    	gxform_args = {'ry',pi/2, 'rz',pi, 'sa',1.25, 'tx',-0.5};
    	    case 'mri',
    		gxform_args = {'rx',pi/2, 'rz',pi, ...
    	    		       'se',size(data_vals) .* 1.25, ...
    	    		       'te',size(data_vals) .* [0.45 0 0.65]};
    	end % switch data

    	   dir_term = {'diry', 0.5, 'dstnt'};
    	  rdir_term = {'diry',-0.1, 'nrby'};

end % switch surf

%%% common terms for SF_DEF_SURF
cterms = {{'batch'}, ...
    	  {'ssize', 'gmean'}, ...
    	  {'datav', 0.5, 1}, ...
    	   dir_term, rdir_term, ...
    	  {'diffe', 0.5}};

%%% display terms for SF_DEF_SURF
if plot_flag
    dterms = {{'surf_f', 5}, {'light', [0 -1 0]}};
else
    dterms = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% common actions start here

disp('computing spectral samples...');

C = sf_gen_intr(fpoints, 'rectq', fpoints, 0.001);
A = sf_gen_area(C,length(C), 'con');
[E,F] = sf_delaunay(C,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% baseline model starts here

disp('beginning baseline model...');

% baseline terms for SF_DEF_SURF (use common and display terms)
bterms = {{'pass', bpass}, cterms{:}, dterms{:}};

%%% generate baseline model and deform
%%% start with rpoints and subsample to get equal ranges
 b0{1}    = sf_gen_intr(rpoints, 'distq', 1);
 b0{2}    = sf_gxform(sf_gen_extr(b0{1}, surf), gxform_args{:});
 b0{3}    = sf_gen_border(b0{1}, 'num_rect');
[b0{4:5}] = sf_delaunay(b0{1}, 4/rpoints);
[b0{1:5}] = sf_sub_surf(b0{1:5}, 'ptot', bpoints.^2);

[bn{1:7}] = sf_def_surf(b0{1:5}, data_vals, bterms);

disp('computing normalized baseline spectrum (for later comparison)...');
Sb = sf_nspec(sf_opt_intr(opt_args{:}, bn{1:5}), C, ...
    	      sf_gen_area(bn{1},1, 'face', bn{5}), A, ...
    	      bn{2}, xform_args);

if xplot_flag
    sf_xplot(C,Sb,F, 'optimized baseline spectra', 'linear', 'mag');
end

clear b0 bn


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% fixed (non-resampling) model starts here

disp('beginning fixed model...');

% non-resampling terms for SF_DEF_SURF (use common and display terms)
nterms = {{'pass', npass}, cterms{:}, dterms{:}, ...
    	  {'errf', 1, C, Sb, 0, xform_args}, ...
    	  };

nerr = [];

for i = 1:reps
    %%% generate non-resampling model and deform
    %%% (start with rpoints and subsample to get equal ranges)
     n0{1}    = sf_gen_intr(rpoints, 'distq', 1);
     n0{2}    = sf_gxform(sf_gen_extr(n0{1}, surf), gxform_args{:});
     n0{3}    = sf_gen_border(n0{1}, 'num_rect');
    [n0{4:5}] = sf_delaunay(n0{1}, 4/rpoints);
    [n0{1:5}] = sf_sub_surf(n0{1:5}, 'ptot', npoints.^2);
    
    [nn{1:7}] = sf_def_surf(n0{1:5}, data_vals, nterms);

    %%% save stats from best pass (minimum spectral error) in nerr
    [val,idx] = min(nn{7}(:,10));
    nerr = [nerr; nn{7}(idx,:)];

    % compute final spectrum (MAY NOT BE SPECTRUM FOR MIN SPEC ERR)
    if xplot_flag
	Sna = sf_nspec(sf_opt_intr(opt_args{:}, nn{1:5}), C, ...
    		       sf_gen_area(nn{1},1, 'face', nn{5}), A, ...
    		       nn{2}, xform_args);
	sf_xplot(C,Sna,F,      'optimized fixed spectra', 'linear', 'mag');
    end

end % for i

clear n0 nn

%%% display results and summarize fixed model
fprintf(fid, ['\nfixed  ', format0]);
fprintf(fid, format1, nerr');

if reps > 1
    nerr_mean      = mean(nerr,1); 
    nerr_mean(2:8) = round(nerr_mean(2:8));
    nerr_std       = std (nerr,1); 
    nerr_std (2:8) = round(nerr_std (2:8));

    fprintf(fid, 'mean of %d reps \n', reps);
    fprintf(fid, format1, nerr_mean);
    fprintf(fid, 'std  of %d reps \n', reps);
    fprintf(fid, format1, nerr_std );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% resampled model starts here

disp('beginning resampled model...');

all_mean = zeros(length(cutofflist), length(threshlist), 2 + length(nerr_mean));
all_std  = zeros(length(cutofflist), length(threshlist), 2 + length(nerr_mean));

%%% test range of possible cutoffs and threshold values
for cutoffind = 1:length(cutofflist)
    for threshind = 1:length(threshlist)
    
	% resampling terms for SF_DEF_SURF (use common and display terms)
	rterms = {{'pass', rpass}, cterms{:}, dterms{:}, ...
    		  {'errf', 1, C, Sb, 0, xform_args}, ...
    		  {'subf', 1, cutofflist(cutoffind), threshlist(threshind), 2, npoints.^2, opt_args, xform_args}, ...
    		  };

	rerr = [];

	for i = 1:reps
	    %%% generate resampling model and deform
	    %%% (start with rpoints and subsample to get equal ranges)
	     r0{1}    = sf_gen_intr(rpoints, 'distq', 1);
	     r0{2}    = sf_gxform(sf_gen_extr(r0{1}, surf), gxform_args{:});
	     r0{3}    = sf_gen_border(r0{1}, 'num_rect');
	    [r0{4:5}] = sf_delaunay(r0{1}, 4/rpoints);

	    [rn{1:7}] = sf_def_surf(r0{1:5}, data_vals, rterms);

	    %%% save stats from best pass (minimum spectral error) in rerr
	    [val,idx] = min(rn{7}(:,10));
	    rerr = [rerr; rn{7}(idx,:)];

	    % compute final spectrum (MAY NOT BE SPECTRUM FOR MIN SPEC ERR)
	    if xplot_flag
		Sra = sf_nspec(sf_opt_intr(opt_args{:}, rn{1:5}), C, ...
    			       sf_gen_area(rn{1},1, 'face', rn{5}), A, ...
    			       rn{2}, xform_args);
		sf_xplot(C,Sra,F,      'optimized resampled spectra', 'linear', 'mag');
	    end

	end % for i

    	%%% display results and summarize resampled model
	fprintf(fid, ['\nresam %s %d %f %f %f %d ', format0], rterms{9}{1:6});
	fprintf(fid, format1, rerr');

	if reps > 1
	    rerr_mean      = mean(rerr,1);
	    rerr_mean(2:8) = round(rerr_mean(2:8));
	    rerr_std       = std (rerr,1);
	    rerr_std (2:8) = round(rerr_std (2:8));

	    fprintf(fid, 'mean of %d reps \n', reps);
	    fprintf(fid, format1, rerr_mean);
	    fprintf(fid, 'std of  %d reps \n', reps);
	    fprintf(fid, format1, rerr_std );

    	    %%% compute statistics relative to fixed model
	    fprintf(fid, 'Prob of significantly different  target  errors (t-test): %4.2f \n', ...
		    sf_ttest(nerr_mean( 9), nerr_std( 9)^2, reps, ...
    			     rerr_mean( 9), rerr_std( 9)^2, reps));
	    fprintf(fid, 'Prob of significantly different spectral errors (t-test): %4.2f \n', ...
		    sf_ttest(nerr_mean(10), nerr_std(10)^2, reps, ...
    			     rerr_mean(10), rerr_std(10)^2, reps));
	end

    	%%% save results to facilitate statistics and plots
    	all_mean(cutoffind, threshind, :) = [ cutofflist(cutoffind) threshlist(threshind) rerr_mean ];
    	all_std (cutoffind, threshind, :) = [ cutofflist(cutoffind) threshlist(threshind) rerr_std  ];
    	
    end % for threshind
end % for cutoffind


fprintf(fid, '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid, '%% %4d/%02d/%02d %02d:%02d:%02d ', fix(clock)); 
fprintf(fid, '\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n');
fclose('all');

more on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

