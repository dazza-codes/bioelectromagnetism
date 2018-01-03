
% test_def_sph.m    - test sf_def_model with spherical data

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES:
% - changed from -depsc (color EPS) to -deps to save disk space

% good parameter combinations:
%   	freq	energy	rpts	npts	npass	rpass
%   	----	------	----	----	-----	-----
%   	0.75    150E-6	10  	20  	15   	10
%   	0.75    150E-6	10  	30  	25   	15
%   	0.8 	200E-6	10  	30  	25  	15
%   	0.75	175E-6	10  	30  	20  	15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% configuration

more off

npass = 20; 	    	% number of passes for non-resampling
rpass = 15; 	    	% number of passes for resampling

fpoints = 30;	    	% number of points for spectral plots
npoints = 30;	    	% number of points for non-resampling
rpoints = 10;	    	% number of points for resampling


  opt_args = {'difeq','deld',0.001};
xform_args = {'traq', 'delf',0.01, 0.1};
 comp_args = {'xyzp', '1d1','1d2','1dm', 'bin','cum', 'csum'};

% generate spherical target surface, restrict to negative x values
[dc,ds,de,df] = sf_gen_surf(40,'isph',1,0.01);
ds = ds(ds(:,1) < 0,:);

cterms =   {{'ssize', 'gmean'}, ...
    	    {'datas', 0.5, 0.05, 1}, ...
    	    {'dirx', 0.5, 'dstnt'}, ...
    	    {'diffe', 0.5}};

dterms =   {{'surf_f', 1}, {'light', [-1 0 0]}};
% dterms = {};

nterms =    {{'pass', npass}, cterms{:}, dterms{:}};
%    	     {'freq_d', 1, xform_args, comp_args}};

rterms =    {{'pass', rpass}, cterms{:}, dterms{:}, ...
    	     {'subf', 1, 0.75, 0.000175, 2, npoints.^2, opt_args, xform_args}, ...
    	     }; 


%%% generate spectral samples
   C  = sf_gen_intr(fpoints, 'rect', fpoints, 0.001);
[E,F] = sf_delaunay(C, 4);
   A  = sf_gen_area(C,length(C), 'con');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deformation without resampling

 n0{1}      = sf_gen_intr(rpoints, 'distq', 1);
 n0{2}	    = sf_gxform(sf_gen_extr(n0{1}, 'sheet'), ...
    	    	        'ry',pi/2, 'rz',pi, 'sa',0.7, 'tx',-0.55);
 n0{3}      = sf_gen_border(n0{1}, 'num_rect');
[n0{4:5}]   = sf_delaunay(n0{1}, 4/rpoints);
[n0{1:5}]   = sf_sub_surf(n0{1:5}, 'ptot', npoints.^2);

[nn{1:7}]   = sf_def_surf(n0{1:5}, ds, nterms);

if (0)
    light('Position', [-1 -1 -1]);
    view([-21 28]);
    caxis([-5 10]); colormap(gray(256)); shading faceted; lighting gouraud;
    % +2 brightness
    
    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.surf-sph-fixed.eps
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
    print -noui -deps fig.apps.spec-sph-fixed.eps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% deformation with resampling

 r0{1}	    = sf_gen_intr(rpoints, 'distq', 1);
 r0{2}	    = sf_gxform(sf_gen_extr(r0{1}, 'sheet'), ...
    	    	    	'ry',pi/2, 'rz',pi, 'sa',0.7, 'tx',-0.55);
 r0{3}	    = sf_gen_border(r0{1}, 'num_rect');
[r0{4:5}]   = sf_delaunay(r0{1}, 4/rpoints);

[rn{1:7}]   = sf_def_surf(r0{1:5}, ds, rterms);


if (0)
    light('Position', [-1 -1 -1]);
    view([-21 28]);
    caxis([-5 10]); colormap(gray(256)); shading faceted; lighting gouraud;
    % +2 brightness

    set(gcf,'PaperPositionMode','auto');
    print -noui -deps fig.apps.surf-sph-resam.eps
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
    print -noui -deps fig.apps.spec-sph-resam.eps
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
