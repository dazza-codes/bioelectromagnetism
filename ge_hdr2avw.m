function [ avw ] = ge_hdr2avw(ge)

% ge_hdr2avw - extract Analyze header from ge struct
% 
% [avw] = ge_hdr2avw(ge)
% 
% ge is a struct obtained from ge_hdr_read
% 
% avw.hdr contains the Analyze header fields.
% 
% There is no reorientation of the volume dimensions 
% in this function, whereas ge_series2avw does 
% reorient the data and corresponding header 
% fields.
% 
% see also ge_series2avw, ge_series_read, ge_hdr_read,
%          avw_hdr_write, avw_hdr_read
%


% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Souheil J. Inati  <souheil.inati@nyu.edu> at 03/2003
% Dartmouth College, May 2000
% 
% Darren.Weber@flinders.edu.au, March 2003
% - Substantially redesigned file handling and function
%   call structures for integration with mri_toolbox at 
%   http://eeg.sf.net
% - Requested permission to distribute code under GPL licence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize the ANALYZE header
avw = avw_hdr_make;

version = '[$Revision: 1.1 $]';
fprintf('\nGE_HDR2AVW [v%s]\n',version(12:16)); tic

% Assign the GE image parameters
avw.hdr.dime.datatype = 4;
avw.hdr.dime.bitpix = ge.hdr.image.screenformat; % usually 16 bits
avw.hdr.dime.glmax = 65535;
avw.hdr.dime.dim(1) =  4;                   % Dimensions (always 4?)
avw.hdr.dime.dim(2) =  ge.hdr.image.imatrix_X; % X Voxels
avw.hdr.dime.dim(3) =  ge.hdr.image.imatrix_Y; % Y Voxels
avw.hdr.dime.dim(4) =  ge.hdr.image.slquant;   % Z Voxels
avw.hdr.dime.dim(5) =  1;                   % Time Points
avw.hdr.dime.pixdim(2) = ge.hdr.image.pixsize_X;                       % X Voxel Size
avw.hdr.dime.pixdim(3) = ge.hdr.image.pixsize_Y;                       % Y Voxel Size
avw.hdr.dime.pixdim(4) = ge.hdr.image.slthick + ge.hdr.image.scanspacing; % Z Voxel Size

t=toc; fprintf('...done (%5.2f sec).\n',t);

return
