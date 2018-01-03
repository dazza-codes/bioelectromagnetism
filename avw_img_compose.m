function avw = avw_img_compose(files,IMGorient,machine),

% AVW_IMG_COMPOSE - Compose single slice Analyze files into a volume
%
% [ avw, machine ] = avw_img_compose(files, orient, machine)
%
% files  - a struct created with F = DIR('*.img').
%          The filenames of the .img files are in F.name and
%          these files are composed into a volume, in the 
%          order given by the DIR command
% 
% orient - force reading IMG in specified orientation, integer values:
%
%          '', read header history orient field
%          0,  transverse/axial unflipped
%          1,  coronal unflipped
%          2,  sagittal unflipped
%          3,  transverse/axial flipped
%          4,  coronal flipped
%          5,  sagittal flipped
%
%          Note that composed volume is given in this orientation
% 
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le' but the routine
%           will automatically switch between little and big
%           endian to read any such Analyze header.  It
%           reports the appropriate machine format and can
%           return the machine value.
%
% Returned values:
% 
% avw.hdr - a struct with image data parameters.
% avw.img - a 3D matrix of image data (double precision).
% 
% See also: AVW_IMG_READ & AVW_HDR_READ (called by this function)
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('files','var'),
    error('AVW_IMG_COMPOSE: No files to compose');
end

if ~exist('IMGorient','var'), IMGorient = ''; end
if ~exist('machine','var'), machine = 'ieee-le'; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of files (slices) to compose
Nfiles = size(files,1);
% Size of first file (slice) in bytes
Fsize = files(1).bytes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise the composed volume

% Read the first slice header
firstslice = avw_hdr_read(files(1).name,machine);

avw = firstslice;

% Not sure if this is always the case, but individual
% slice files output by http://xmedcon.sourceforge.net
% have the FOV in .hdr.dime.pixdim(4) after reading
% some Siemens .ima files.  It might be worth checking
% at this point and dividing this FOV by the
% number of slices
if avw.hdr.dime.pixdim(4) > 20,
    msg = sprintf('AVW_IMG_COMPOSE: slice pixdim(4) is very large, assuming it is FOV and converting!\n');
    warning(msg);
    % OK this field is probably FOV, so lets
    % divide it by the total number of slices
    SliceFOV = double(avw.hdr.dime.pixdim(4));
    avw.hdr.dime.pixdim(4) = single(SliceFOV / Nfiles);
end

% Initialise avw.img and reset some header dimensions
avw.hdr.dime.dim(4) = Nfiles;
avw = avw_img_init(avw);

% Now reset orient field of avw, as all data will
% be in standard axial unflipped orientation after
% avw_img_read.
avw.hdr.hist.orient = char(0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add slices, depends on whether slice order is normal or flipped

orient = double(firstslice.hdr.hist.orient);

if orient < 3,
    % process slices as unflipped volume
    for f = 1:Nfiles,
        if files(f).bytes == Fsize,
            slice = avw_img_read(files(f).name,IMGorient,machine);
            avw   = avw_add_slice(avw,slice,f,orient);
        else
            msg = sprint('AVW_IMG_COMPOSE: This file is not the same size as the first file!');
            error(msg);
        end
    end
else
    % process slices as flipped volume
    for f = Nfiles:-1:1,
        if files(f).bytes == Fsize,
            slice = avw_img_read(files(f).name,IMGorient,machine);
            avw   = avw_add_slice(avw,slice,f,orient);
        else
            msg = sprint('AVW_IMG_COMPOSE: This file is not the same size as the first file!');
            error(msg);
        end
    end
end


return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avw] = avw_add_slice(avw,slice,Nslice,orient),
    
    % Check that composed volume and current slice
    % have the same orientation
    Vorient = double(  avw.hdr.hist.orient);
    Sorient = double(slice.hdr.hist.orient);
    
    if Vorient ~= Sorient,
        msg = sprintf('AVW_IMG_COMPOSE: Slice %3d has different orientation to volume',Nslice);
        error(msg);
    end
    
    switch double(orient),
    
    case 0, % axial unflipped
        
        % Slices in 'z' axis - from patient inferior to superior
        avw.img(:,:,Nslice) = slice.img;
        
    case 1, % coronal unflipped
        
        % Slices in 'y' axis - from patient posterior to anterior
        avw.img(:,Nslice,:) = slice.img;
        
    case 2, % sagittal unflipped
        
        % Slices in 'x' axis - from patient right to left
        avw.img(Nslice,:,:) = slice.img;
        
    case 3, % axial flipped
        
        % Slices in 'z' axis - from patient inferior to superior
        avw.img(:,:,Nslice) = slice.img;
        
    case 4, % coronal flipped
        
        % Slices in 'y' axis - from patient anterior to posterior
        avw.img(:,Nslice,:) = slice.img;
        
    case 5, % sagittal flipped
        
        % Slices in 'x' axis - from patient right to left
        avw.img(Nslice,:,:) = slice.img;
    end
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [avw] = avw_img_init(avw),
    
	PixelDim = double(avw.hdr.dime.dim(2));
	RowDim   = double(avw.hdr.dime.dim(3));
	SliceDim = double(avw.hdr.dime.dim(4));
    
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    
    switch double(avw.hdr.hist.orient),
    
    case {0,3}, % axial unflipped or flipped
        
        x = PixelDim;   Xsz = PixelSz;
        y = RowDim;     Ysz = RowSz;
        z = SliceDim;   Zsz = SliceSz;
        
    case {1,4}, % coronal unflipped or flipped
        
        x = PixelDim;   Xsz = PixelSz;
        y = SliceDim;   Ysz = SliceSz;
        z = RowDim;     Zsz = RowSz;
        
    case {2,5}, % sagittal unflipped or flipped
        
        x = SliceDim;   Xsz = SliceSz;
        y = PixelDim;   Ysz = PixelSz;
        z = RowDim;     Zsz = RowSz;
        
    otherwise, % assume axial unflipped or flipped
        
        msg = sprintf('AVW_IMG_COMPOSE: No specified orientation\n');
        warning(msg);
        
        x = PixelDim;   Xsz = PixelSz;
        y = RowDim;     Ysz = RowSz;
        z = SliceDim;   Zsz = SliceSz;
        
    end
    
    avw.img = zeros(x,y,z);
    avw.hdr.dime.dim(2) = x;
	avw.hdr.dime.dim(3) = y;
	avw.hdr.dime.dim(4) = z;
    
    avw.hdr.dime.pixdim(2) = Xsz;
    avw.hdr.dime.pixdim(3) = Ysz;
    avw.hdr.dime.pixdim(4) = Zsz;
    
return
