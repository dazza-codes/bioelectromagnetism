function avw_write(avw, fileprefix, IMGorient, machine, verbose)

% avw_write - write Analyze files (*.img & *.hdr)
% 
% This function is a wrapper for
% avw_hdr_write and avw_img_write.
% 
% avw_write(avw,fileprefix,[IMGorient],[machine],[verbose])
% 
% where the input avw is a struct (see avw_read):
% 
% avw.hdr    - a struct with image data parameters.
% avw.img    - a 3D matrix of image data (double precision).
% 
% fileprefix - a string, the filename without the .img
%              extension. If empty, may use avw.fileprefix
% 
% IMGorient - optional int; force writing of specified 
%             orientation, as follows:
% 
%   [],  if empty, will use avw.hdr.hist.orient field
%    0,  transverse/axial unflipped (default, radiological)
%    1,  coronal unflipped
%    2,  sagittal unflipped
%    3,  transverse/axial flipped
%    4,  coronal flipped
%    5,  sagittal flipped
% 
% This function will set avw.hdr.hist.orient and write the 
% image data in a corresponding order.  This function is in
% alpha development, it has not been exhaustively tested 
% (as of 07/2003).
% 
% See the comments and notes throughout avw_hdr_read and 
% avw_img_read for more information and documentation on 
% the Analyze orientation options.  Orientations 3-5 are
% NOT recommended! They are part of the Analyze format, 
% but only used in Analyze for faster raster graphics 
% during movies. (Also see the copy of the Analyze 7.5 
% format pdf in the mri_toolbox doc folder; although 
% the code is largely self sufficient, that pdf is the 
% authoritative description).
% 
% machine - optional string; see machineformat in fread
%           for details. The default here is 'ieee-le'.
% 
% verbose - the default is to output processing information to
%           the command window.  If verbose = 0, this will not 
%           happen.
%
% Tip: to change the data type, set avw.hdr.dime.datatype to:
% 
%     1    Binary             (  1 bit  per voxel)
%     2    Unsigned character (  8 bits per voxel)
%     4    Signed short       ( 16 bits per voxel)
%     8    Signed integer     ( 32 bits per voxel)
%    16    Floating point     ( 32 bits per voxel)
%    32    Complex, 2 floats  ( 64 bits per voxel), not supported
%    64    Double precision   ( 64 bits per voxel)
%   128    Red-Green-Blue     (128 bits per voxel), not supported
% 
% See also: avw_read, avw_hdr_write, avw_img_write
%           avw_hdr_read, avw_img_read 
%           avw_view
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze 7.5 format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%                    This code attempts to respect the integrity of
%                    the format, although no guarantee is given that
%                    it has been exhaustively tested for compatibility
%                    with Analyze software (it should be though).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var'),
    doc avw_write;
    msg = sprintf('\n...no input avw.\n');
    error(msg);
elseif isempty(avw),
    msg = sprintf('\n...empty input avw.\n');
    error(msg);
elseif ~isfield(avw,'img'),
    msg = sprintf('\n...empty input avw.img\n');
    error(msg);
end

if ~exist('fileprefix','var'),
    if isfield(avw,'fileprefix'),
        if ~isempty(avw.fileprefix),
            fileprefix = avw.fileprefix;
        else
            fileprefix = [];
        end
    else
        fileprefix = [];
    end
end
if isempty(fileprefix),
    [fileprefix, pathname, filterindex] = uiputfile('*.hdr','Specify an output Analyze .hdr file');
    if pathname, cd(pathname); end
    if ~fileprefix,
        doc avw_write;
        error('no output .hdr file specified');
    end
end

if findstr('.hdr',fileprefix),
%    fprintf('AVW_WRITE: Removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if findstr('.img',fileprefix),
%    fprintf('AVW_WRITE: Removing .img extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.img','');
end

if ~exist('IMGorient','var'), IMGorient = ''; end
if ~exist('machine','var'), machine = 'ieee-le'; end
if ~exist('verbose','var'), verbose = 1; end

if isempty(IMGorient), IMGorient = ''; end
if isempty(machine), machine = 'ieee-le'; end
if isempty(verbose), verbose = 1; end


%------------------------------------------------------------------------
% MAIN

version = '[$Revision: 1.1 $]';

avw_img_write(avw,fileprefix,IMGorient,machine,verbose);

% MUST write header after the image, to ensure any
% orientation changes during image write are saved
% in the header.  The header is saved by avw_img_write, so it's
% not necessary to call it here.  If it were called here, it
% would be called as such:
% avw_hdr_write(avw,fileprefix,machine,verbose);

return
