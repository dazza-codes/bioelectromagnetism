
% test_def_mri.m    - test sf_def_model with MRI data

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES:
% - changed from -depsc (color EPS) to -deps to save disk space
    	    	    	
% good parameter combinations:
%   	freq	energy	rpts	npts	npass	rpass
%   	----	------	----	---- 	-----	-----
%   	0.75    150E-6	10  	20  	20   	15
%   	0.75    150E-6	10  	30  	35   	20
%   	0.75	175E-6	10  	30  	35  	25
%   	0.80	175E-6	10  	30  	30    	20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% configuration

more off

npass = 30; 	    % number of passes for non-resampling
rpass = 20; 	    % number of passes for resampling

fpoints = 30;	    % number of points for spectral plots
npoints = 30;	    % number of points for non-resampling
rpoints = 10;	    % number of points for resampling


  opt_args = {'difcq','deld',0.001};	    % arguments for sf_opt_intr
xform_args = {'traq', 'delf',0.01, 0.1};    % arguments for sf_xform
 comp_args = {'xyzp','2d'}; 	    	    % arguments for sf_comp

%%% read data vals
load clkmri003
data_vals = sf_find_3d(clkmri003,80,'vol');
clear clkmri003

gxform_args = {'rx',pi/2, 'rz',pi, ...
    	       'se',size(data_vals) .* 1.25, ...
    	       'te',size(data_vals) .* [0.45 0 0.65]};

cterms = {{'ssize', 'gmean'}, ...
    	  {'datav', 0.5, 1}, ...
    	  {'diry', 0.5, 'dstnt'}, ...
    	  {'diry',-0.05, 'nrby' }, ...
    	  {'diffe', 0.5}};

dterms = {{'surf_f', 1}, {'light', [0 -1 0]}};
% dterms = {};

nterms = {{'pass', npass}, cterms{:}, dterms{:}};

rterms = {{'pass', rpass}, cterms{:}, dterms{:}, ...
    	  {'subf', 1, 0.80, 0.000175, 2, npoints.^2, opt_args, xform_args}, ...
    	  };

%%% generate spectral samples
   C  = sf_gen_intr(fpoints, 'rectq', fpoints, 0.001);
[E,F] = sf_delaunay(C,4);
   A  = sf_gen_area(C,length(C), 'con');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deformation without resampling

 n0{1}    = sf_gen_intr(rpoints, 'distq', 1);
 n0{2}    = sf_gxform(sf_gen_extr(n0{1}, 'sheet'), gxform_args{:});
 n0{3}    = sf_gen_border(n0{1}, 'num_rect');
[n0{4:5}] = sf_delaunay(n0{1}, 4/rpoints);
[n0{1:5}] = sf_sub_surf(n0{1:5}, 'ptot', npoints.^2);

[nn{1:7}] = sf_def_surf(n0{1:5}, data_vals, nterms);

if (0)
    light('Position', [-1 -1 -1]);
    view([-21 28]);
    caxis([-5 10]); colormap(gray(256)); shading faceted; lighting gouraud;
    % +2 brightness

    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.surf-mri-fixed.eps
end


Sna = sf_nspec(sf_opt_intr(opt_args{:}, nn{1:5}), C, ...
    	       sf_gen_area(nn{1},1, 'face', nn{5}), A, ...
    	       nn{2}, xform_args);
sf_xplot(C,Sna,F, 'XYZ', 'linear', 'mag');

if (0)
    % set colormap to gray256, caxis to [0.0 0.1], 
    % axis labels to f1f2, range to -15 15
    % +4 brightness

    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.spec-mri-fixed.eps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deformation with resampling

 r0{1}    = sf_gen_intr(rpoints, 'distq', 1);
 r0{2}    = sf_gxform(sf_gen_extr(r0{1}, 'sheet'), gxform_args{:});
 r0{3}    = sf_gen_border(r0{1}, 'num_rect');
[r0{4:5}] = sf_delaunay(r0{1}, 4/rpoints);

[rn{1:7}] = sf_def_surf(r0{1:5}, data_vals, rterms);


if (0)
    light('Position', [-1 -1 -1]);
    view([-21 28]);
    caxis([-5 10]); colormap(gray(256)); shading faceted; lighting gouraud;
    % +2 brightness

    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.surf-mri-resam.eps
end


Sra = sf_nspec(sf_opt_intr(opt_args{:}, rn{1:5}), C, ...
    	       sf_gen_area(rn{1},1, 'face', rn{5}), A, ...
    	       rn{2}, xform_args);
sf_xplot(C,Sra,F, 'XYZ', 'linear', 'mag');

if (0)
    % set colormap to gray256, caxis to [0.0 0.1],
    % axis labels to f1f2, range to -15 15
    % +4 brightness

    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.spec-mri-resam.eps
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks


