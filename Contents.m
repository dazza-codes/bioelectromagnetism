% MRI_TOOLBOX v2.0
%
% Copyright (C) 2004  Darren Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%  avw_read - read Analyze format data image (*.img)
%  avw_hdr_read - read Analyze format data header (*.hdr)
%  avw_img_read - read Analyze format data image (*.img)
%  avw_img_compose - compose Analyze volume from single slices
%
%  avw_hdr_check_datatype - read Analyze format data header (*.hdr)
%  avw_hdr_make - create Analyze format data header (avw.hdr)
%
%  avw_write - write Analyze files (*.img & *.hdr)
%  avw_hdr_write - write Analyze header file (*.hdr)
%  avw_img_write - write Analyze image files (*.img)
%
%  avw_view - create and navigate ortho views of Analyze 7.5 volume
%  avw_view_hdr - view and modify Analyze header file
%
%  avw_roi - extract a region of interest from avw.img
%  avw_stats - calculate region of interest stats for avw data
%  avw_metric - convert image location from CRS to meters
%
%  avw_binary - return 0 for avw.img <= thresh, 1 otherwise
%  avw_smooth - Guassian smoothing
%
%  avw_center - find center of a volume
%  avw_flip - ortho-flip Analyze data image (avw.img)
%
%  avw2brainstorm - Convert Analyze struct into BrainStorm file
%  avw2cor - converts an avw struct to FreeSurfer COR-* files
%  cor2avw - Read Freesurfer MRI data (COR-001 to COR-256)
%  cor_img_read - Read Freesurfer format data (COR-001 to COR-256)
%
%  mni2tal - MNI to Talairach coordinates (best guess)
%  mni2tal_matrix - Talairach to MNI coordinates (best guess)
%  tal2mni - Talairach to MNI coordinates
%
%  avw_converter_script - convert byte order of .img files
%  avw_homogeneity_correction - 2D correction of MRI RF nonuniformity
%
%  emse_elec2mri - convert point in head frame to mri coordinates
%  emse_mri2elec - convert mri coordinates to points in head frame
%  emse_open_reg - read EMSE/MRVU coregistration matrices
%
%  ge_hdr2avw - extract Analyze header from ge struct
%  ge_hdr_read - read the header info from a GE LX2 or 5.X file
%  ge_series2avw - convert a GE series to Analyze
%  ge_series_read - read a volume of images from a GE series
%
%  ctf_write_mri - write a CTF .mri file
%  ctf_make_mri - create a CTF mri struct with zero image
%  ctf_read_mri - read a CTF .mri file
%
%  gui_avw_open - Load & Display Analyze 7.5 data
%  gui_cor_open - Load & Display FreeSurfer COR-??? data
%  gui_ge_open - Load & Display MRI data
%  gui_mri_open - Load & Display MRI data
%  mri_open - function to call various mri data tools
%  mri_toolbox_recent - Keep track of mri_toolbox .mat files
%  mri_updateparent - General GUI data handing for MRI Toolbox
%
