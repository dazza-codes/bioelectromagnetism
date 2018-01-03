% Spatial Frequency Toolbox
% started 6/96 by Clif Kussmaul
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% sf_comp		- create composite surface spectrum estimates
% sf_def_surf		- deform surface to match target data
% sf_delaunay		- create edges and faces using Delaunay triangulation
% sf_filter		- filter extrinsic values using intrinsic values
% sf_find_3d		- find specific values in 3D data set
% sf_find_neighbors	- find neighbors of given points in target points
% sf_ftest		- F-TEST for significantly different variances
% sf_gen_area		- create area estimate
% sf_gen_border		- identify border points from coordinate values
% sf_gen_edge		- create edges (endpoint indices) from points
% sf_gen_extr		- create extrinsic surface points
% sf_gen_face		- create faces (corner indices) from edges
% sf_gen_intr 		- create intrinsic surface points
% sf_gen_map		- create mapping functions between representations
% sf_gen_surf		- create intr and extr points for special cases
% sf_gxform		- apply geometric transforms (rotate,scale,translate)
% sf_hist		- plot histograms of surface values
% sf_norm		- compute per-column norms
% sf_opt_intr		- optimize intrinsic points
% sf_plot		- plot surface (points/edges/faces) (movie)
% sf_profile		- profile command with specified arguments
% sf_read_3d		- read 3D imaging data from file(s)
% sf_read_byu		- read MOVIE.BYU data from file
% sf_sample		- return random sampling of values
% sf_save_mpg		- save movie as MPG file
% sf_save_vrml		- save surface as VRML file
% sf_show		- show orthogonal slices from 3D data
% sf_sphharm		- create/plot spherical harmonics
% sf_splot		- create movie of surface grid (interp from points)
% sf_spy		- enhanced SPY for full and sparse matrices
% sf_stat		- compute pointwise statistics across signals
% sf_sub_surf	    	- subdivide surface (new points, edges, faces)
% sf_tool		- supporting tools (mostly user interface)
% sf_ttest		- T-TEST for significantly different means
% sf_xform		- spatial frequency transform (multiple methods)
% sf_xplot		- create cross plots for extrinsic points
%

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (display help message if this m-file is executed)
help Contents
