function avw_img_write(avw, fileprefix, IMGorient, machine)

% avw_img_write - write Analyze image files (*.img)
% 
% avw_img_write(avw, fileprefix, [IMGorient], [machine])
% 
% avw.img    - a 3D matrix of image data (double precision).
% avw.hdr    - a struct with image data parameters.  If
%              not empty, this function calls avw_hdr_write.
% 
% fileprefix - a string, the filename without the .img
%              extension. If empty, may use avw.fileprefix
% 
% IMGorient - optional int, force writing of specified 
%             orientation, with values:
% 
%   [],  if empty, will use avw.hdr.hist.orient field
%    0,  transverse/axial unflipped (default, radiological)
%    1,  coronal unflipped
%    2,  sagittal unflipped
%    3,  transverse/axial flipped, left to right
%    4,  coronal flipped, anterior to posterior
%    5,  sagittal flipped, superior to inferior
% 
% This function will set avw.hdr.hist.orient and write the 
% image data in a corresponding order.  This function is 
% in alpha development, so it has not been exhaustively 
% tested (07/2003). See avw_img_read for more information 
% and documentation on the orientation option.  
% Orientations 3-5 are NOT recommended!  They are part 
% of the Analyze format, but only used in Analyze
% for faster raster graphics during movies.
% 
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le'.
% 
% Tip: to change the data type, set avw.hdr.dime.datatype to:
% 
%     1    Binary             (  1 bit  per voxel)
%     2    Unsigned character (  8 bits per voxel)
%     4    Signed short       ( 16 bits per voxel)
%     8    Signed integer     ( 32 bits per voxel)
%    16    Floating point     ( 32 bits per voxel)
%    32    Complex, 2 floats  ( 64 bits per voxel), flipping not supported
%    64    Double precision   ( 64 bits per voxel)
%   128    Red-Green-Blue     (128 bits per voxel), not supported
% 
% See also: avw_write, avw_hdr_write, 
%           avw_read, avw_hdr_read, avw_img_read, avw_view
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2002, Darren.Weber@flinders.edu.au
%                    The Analyze format is copyright 
%                    (c) Copyright, 1986-1995
%                    Biomedical Imaging Resource, Mayo Foundation
%           07/2004, chodkowski@kennedykrieger.org, added ability to
%                    write complex AVW .img files.  added error if complex
%  data is to be flipped.  currently there is no logic for flipping
%  complex data.  also force invalid IMGorient to 0, ie, no flipping.
%  forcing an invalid IMGorient to 0 allows me to remove duplicate code
%  in the IMGorient case/otherwise logic.  i also pulled the fwrite of
%  non-flipped data, IMGorient == 0, out of any looping mechanism.  looping 
%  is not necessary as the data is already in its correct orientation.  
%  removing fwrite from the loops should be faster but, more importantly, 
%  it allows the writing of complex data which contains two matrix 
%  elements per one pixel.  moved the fclose statement from write_image 
%  to avw_img_write, the same function that calls fopen.
%
%  write complex data example:
%     % assume rr contains n-dimensional real values
%     % assume ii contains n-dimensional imaginary values
%     % where the magnitude image can be computed as follows:
%
%     magn = sqrt( rr.^2 + ii.^2 );
%
%     cc = zeros( prod( size( rr ))*2, 1 );     % "*2" bc 1pix = [real, imag]
%     cc(1:2:end) = reshape( rr, prod(size(rr)), 1 );
%     cc(2:2:end) = reshape( ii, prod(size(ii)), 1 );
%
%     avw = avw_hdr_make;
%
%     % numDims, xdim, ydim, zdim, ..., padded with 8 0s
%     tmp = [ ndims(rr) size(rr) zeros(1,8) ];
%     avw.hdr.dime.dim = tmp(1:8);         % keep only the first 8 values
%
%     avw.hdr.dime.datatype = 32;          % complex
%     avw.hdr.dime.bitpix   = 64;          % 4+4 bytes/pixel
%
%     avw.img = cc;
%     avw_img_write( avw, 'foo' );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------
% Check inputs

if ~exist('avw','var'),
    doc avw_img_write;
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
            fileprefix = '';
        end
    end
end
if isempty(fileprefix),
    doc avw_img_write;
    fprintf('\n...no input fileprefix - see help avw_img_write\n');
    return;
end

if ~exist('IMGorient','var'), IMGorient = ''; end
if ~exist('machine','var'), machine = 'ieee-le'; end

if findstr('.hdr',fileprefix),
%    fprintf('AVW_IMG_WRITE: Removing .hdr extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.hdr','');
end
if findstr('.img',fileprefix),
%    fprintf('AVW_IMG_WRITE: Removing .img extension from ''%s''\n',fileprefix);
    fileprefix = strrep(fileprefix,'.img','');
end



%------------------------------------------------------------------------
% MAIN

version = '[$Revision: 1.1 $]';
fprintf('\nAVW_IMG_WRITE [v%s]\n',version(12:16));  tic;

fid = fopen(sprintf('%s.img',fileprefix),'w',machine);
if fid < 0,
    msg = sprintf('Cannot open file %s.img\n',fileprefix);
    error(msg);
else
    avw = write_image(fid,avw,fileprefix,IMGorient,machine);
end
fclose(fid);

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

% MUST write header after the image, to ensure any
% orientation changes during image write are saved
% in the header
avw_hdr_write(avw,fileprefix,machine);

return




%-----------------------------------------------------------------------------------

function avw = write_image(fid,avw,fileprefix,IMGorient,machine)

% short int bitpix;    /* Number of bits per pixel; 1, 8, 16, 32, or 64. */ 
% short int datatype   /* Datatype for this image set */
% /*Acceptable values for datatype are*/ 
% #define DT_NONE             0
% #define DT_UNKNOWN          0    /*Unknown data type*/ 
% #define DT_BINARY           1    /*Binary             ( 1 bit per voxel)*/ 
% #define DT_UNSIGNED_CHAR    2    /*Unsigned character ( 8 bits per voxel)*/ 
% #define DT_SIGNED_SHORT     4    /*Signed short       (16 bits per voxel)*/ 
% #define DT_SIGNED_INT       8    /*Signed integer     (32 bits per voxel)*/ 
% #define DT_FLOAT           16    /*Floating point     (32 bits per voxel)*/ 
% #define DT_COMPLEX         32    /*Complex,2 floats   (64 bits per voxel)/* 
% #define DT_DOUBLE          64    /*Double precision   (64 bits per voxel)*/ 
% #define DT_RGB            128    /*A Red-Green-Blue datatype*/
% #define DT_ALL            255    /*Undocumented*/

switch double(avw.hdr.dime.datatype),
case   1,
    avw.hdr.dime.bitpix = int16( 1); precision = 'bit1';
case   2,
    avw.hdr.dime.bitpix = int16( 8); precision = 'uchar';
case   4,
    avw.hdr.dime.bitpix = int16(16); precision = 'int16';
case   8,
    avw.hdr.dime.bitpix = int16(32); precision = 'int32';
case  16,
    avw.hdr.dime.bitpix = int16(32); precision = 'single';
case  32,
    avw.hdr.dime.bitpix = int16(64); precision = 'float';
case  64,
    avw.hdr.dime.bitpix = int16(64); precision = 'double';
case 128,
    error('...RGB datatype not yet supported.\n');
otherwise
    fprintf('...unknown datatype, using type 16 (32 bit floats).\n');
    avw.hdr.dime.datatype = int16(16);
    avw.hdr.dime.bitpix = int16(32); precision = 'single';
end


% write the .img file, depending on the .img orientation
fprintf('...writing %s precision Analyze image (%s).\n',precision,machine);

fseek(fid,0,'bof');

% The standard image orientation is axial unflipped
if isempty(avw.hdr.hist.orient),
    msg = [ '...WARNING: avw.hdr.hist.orient ~= 0.\n',...
            '   This function assumes the input avw.img is\n',...
            '   in axial unflipped orientation in memory.  This is\n',...
            '   created by the avw_img_read function, which converts\n',...
            '   any input file image to axial unflipped in memory.\n'];
    fprintf(msg)
end

if isempty(IMGorient),
    fprintf('...no IMGorient specified, using avw.hdr.hist.orient value.\n');
    IMGorient = double(avw.hdr.hist.orient);
end

if ~isfinite(IMGorient),
    fprintf('...WARNING: IMGorient is not finite!\n');
    fprintf('...unknown orientation specified, assuming default axial unflipped\n');
    fprintf('...using case 0 \n' );
    IMGorient = 0;
end

maxCase = 5;
if (IMGorient > maxCase ),
    fprintf('...WARNING: IMGorient is greater than %d!\n', maxCase);
    fprintf('...unknown orientation specified, assuming default axial unflipped\n');
    fprintf('...using case 0 \n' );
    IMGorient = 0;
end

% if datatype is complex (32), do not allow flipping bc IMGorient == 0 is
% the only fwrite logic that will handle 2 elements/pixel.  i modified the
% non-flipping fwrite logic to write out then entire avw.img in a single
% block, ie, fwrite does not care about the volume dimensions.  this allows
% one to write out complex data where one pixel is represented by 2 elements
% in the avw.img matrix.  note that AVW complex data means:
%     [ [r1,i1], [r2,i2], ... ]

% issue an error if complex data needs to be flipped
if (( avw.hdr.dime.datatype == 32 ) && ( IMGorient ~= 0 ))
   msg = [ '...ERROR:  avw.hdr.dime.datatype = 32 (complex) and IMGorient, ', ...
           'the orientation of the volume, is not 0.  A non-zero ', ...
           'IMGorient requires flipping the data.  Flipping is not ', ...
           'implemented for complex data.  Flip your data before ', ...
           'calling this function' ];
   msg = sprintf( '%s (%s).', msg, mfilename );
   error( msg );
end;


switch IMGorient,
    
case 0, % transverse/axial unflipped
    
    % For the 'transverse unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient posterior to anterior
    % Slices in 'z' axis                  - from patient inferior to superior
    
    fprintf('...writing axial unflipped\n');
    
    avw.hdr.hist.orient = uint8(0);
    
    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    newWay = true;
    if ( newWay )
       % since there is no flipping, write out the entire volume at once
       fwrite( fid, avw.img, precision );
    else
       x = 1:PixelDim;
       for z = 1:SliceDim,
           for y = 1:RowDim,
               fwrite(fid,avw.img(x,y,z),precision);
           end
       end
    end;
    
case 1, % coronal unflipped
    
    % For the 'coronal unflipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'y' axis                  - from patient posterior to anterior
    
    fprintf('...writing coronal unflipped\n');
    
    avw.hdr.hist.orient = uint8(1);
    
    SliceDim = double(avw.hdr.dime.dim(3)); % y
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(3));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for y = 1:SliceDim,
        for z = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 2, % sagittal unflipped
    
    % For the 'sagittal unflipped' type, the voxels are stored with
    % Pixels in 'y' axis (varies fastest) - from patient posterior to anterior
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'x' axis                  - from patient right to left
    
    fprintf('...writing sagittal unflipped\n');
    
    avw.hdr.hist.orient = uint8(2);
    
    SliceDim = double(avw.hdr.dime.dim(2)); % x
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(3)); % y
    SliceSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(3));
    
    y = 1:PixelDim;
    for x = 1:SliceDim,
        for z = 1:RowDim,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 3, % transverse/axial flipped
    
    % For the 'transverse flipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'y' axis                  - from patient anterior to posterior*
    % Slices in 'z' axis                  - from patient inferior to superior
    
    fprintf('...writing axial flipped (+Y from Anterior to Posterior)\n');
    
    avw.hdr.hist.orient = uint8(3);
    
    SliceDim = double(avw.hdr.dime.dim(4)); % z
    RowDim   = double(avw.hdr.dime.dim(3)); % y
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(4));
    RowSz    = double(avw.hdr.dime.pixdim(3));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for z = 1:SliceDim,
        for y = RowDim:-1:1, % flipped in Y
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 4, % coronal flipped
    
    % For the 'coronal flipped' type, the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient inferior to superior
    % Slices in 'y' axis                  - from patient anterior to posterior
    
    fprintf('...writing coronal flipped (+Z from Superior to Inferior)\n');
    
    avw.hdr.hist.orient = uint8(4);
    
    SliceDim = double(avw.hdr.dime.dim(3)); % y
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(2)); % x
    SliceSz  = double(avw.hdr.dime.pixdim(3));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(2));
    
    x = 1:PixelDim;
    for y = 1:SliceDim,
        for z = RowDim:-1:1,
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    
case 5, % sagittal flipped
    
    % For the 'sagittal flipped' type, the voxels are stored with
    % Pixels in 'y' axis (varies fastest) - from patient posterior to anterior
    % Rows in   'z' axis                  - from patient superior to inferior
    % Slices in 'x' axis                  - from patient right to left
    
    fprintf('...writing sagittal flipped (+Z from Superior to Inferior)\n');
    
    avw.hdr.hist.orient = uint8(5);
    
    SliceDim = double(avw.hdr.dime.dim(2)); % x
    RowDim   = double(avw.hdr.dime.dim(4)); % z
    PixelDim = double(avw.hdr.dime.dim(3)); % y
    SliceSz  = double(avw.hdr.dime.pixdim(2));
    RowSz    = double(avw.hdr.dime.pixdim(4));
    PixelSz  = double(avw.hdr.dime.pixdim(3));
    
    y = 1:PixelDim;
    for x = 1:SliceDim,
        for z = RowDim:-1:1, % superior to inferior
            fwrite(fid,avw.img(x,y,z),precision);
        end
    end
    

%% this fall-thru case should never happen since an unknown IMGorient
%% is detected before this case statement and forced to "0".  no need
%% to have duplicate logic as this block of code is the same as case 0.

%otherwise, % transverse/axial unflipped
%    
%    % For the 'transverse unflipped' type, the voxels are stored with
%    % Pixels in 'x' axis (varies fastest) - from patient right to left
%    % Rows in   'y' axis                  - from patient posterior to anterior
%    % Slices in 'z' axis                  - from patient inferior to superior
%    
%    fprintf('...unknown orientation specified, assuming default axial unflipped\n');
%    
%    avw.hdr.hist.orient = uint8(0);
%    
%    SliceDim = double(avw.hdr.dime.dim(4)); % z
%    RowDim   = double(avw.hdr.dime.dim(3)); % y
%    PixelDim = double(avw.hdr.dime.dim(2)); % x
%    SliceSz  = double(avw.hdr.dime.pixdim(4));
%    RowSz    = double(avw.hdr.dime.pixdim(3));
%    PixelSz  = double(avw.hdr.dime.pixdim(2));
%    
%    x = 1:PixelDim;
%    for z = 1:SliceDim,
%        for y = 1:RowDim,
%            fwrite(fid,avw.img(x,y,z),precision);
%        end
%    end
    
end

% Update the header
avw.hdr.dime.dim(2:4) = int16([PixelDim,RowDim,SliceDim]);
avw.hdr.dime.pixdim(2:4) = single([PixelSz,RowSz,SliceSz]);

return
