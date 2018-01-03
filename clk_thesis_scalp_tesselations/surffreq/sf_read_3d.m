function data = sf_read_3d(filename, scale)
% SF_READ_3D	Read 3D imaging data from file(s)
%   	    	data = sf_read_3d(filename [, scale])
%   	    	using SCALE to adjust amplitudes

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - to look at one slice
%   image(squeeze(data(:,:,30)))

% - to look at collapsed slices

%%% THINGS TO DO
% ? duplicate or alternate data values to get square voxels
% ? option to display slices as they are read
% ? more formats (*.anat 1147818, 12897362, 5953141)
% ? fill structure with useful values from file headers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments 
if (nargin < 1)
    help sf_read_3d; return;
elseif ~ischar(filename)
    error(['FILENAME must be a string.']);
elseif (nargin < 2)
    scale = 1;
end


%%% read from file depending on size
[fid size] = sf_read_3d_open(filename, 'r');

if     size < 131072
    %%% list of images
    fcnt = 0;
    while 1
    	fname = fgetl(fid);
    	if ~ischar(fname); break; end
    	fcnt = fcnt + 1;
    	flist{fcnt} = fname;
    end
    data = [];
    for j=1:fcnt
    	data = cat(3, data, sf_read_3d(flist{j}, scale));
    end

elseif size == 138976
    data = sf_read_3d_read(fid, 7904, [256 256], 'int16') * scale / 256;
    data = sf_read_3d_truncate(data);
    
elseif size == 532192
    data = sf_read_3d_read(fid, 7904, [512 512], 'int16') * scale / 256;
    data = sf_read_3d_truncate(data);    

elseif size == 145408
    %%% 256x256 GE MRI file
    data = sf_read_3d_read(fid, 14336, [256 256], 'int16') * scale;
    data = sf_read_3d_truncate(data);

elseif size == 538624
    %%% 512x512 GE MRI file
    data = sf_read_3d_read(fid, 14336, [512 512], 'int16') * scale;
    data = sf_read_3d_truncate(data);

elseif size == 135168
    %%% 256x256 Siemens MRI file (little endian)
    fclose(fid); fid = fopen(filename, 'r', 'l');
    data = sf_read_3d_read(fid,  4096, [256 256], 'int16') * scale / 256;
    data = sf_read_3d_truncate(data);
    
elseif size == 528384
    %%% 512x512 Siemens MRI file (little endian)
    fclose(fid); fid = fopen(filename, 'r', 'l');
    data = sf_read_3d_read(fid,  4096, [512 512], 'int16') * scale / 256;
    data = sf_read_3d_truncate(data);

elseif size == 1565696
    %%% Siemens PET file 
    fclose(fid); fopen(filename, 'r', 'l');
    data = [];
    for j = 1:46
    	offset = 2048 + j * (512 + 128*128*2);
    	data1 = sf_read_3d_read(fid, offset, [128 128], 'int16') * scale / 256;
   	data = cat(3, data, sf_read_3d_truncate(data1));
    end
    
else
    error(['Unknown file size: ', num2str(size), '.']);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fid, size] = sf_read_3d_open(filename, mode)
% SF_READ_3D_OPEN   	open file and return fid and size

disp(['Reading ',filename, '...']);
% try to open file
fid = fopen(filename, mode);
if fid < 0; error(['Unable to open ', filename, '.']); end

if nargout > 1
    % seek to EOF to determine size
    status = fseek(fid, 0, 'eof');
    if status == -1; error(ferror(fid)); end
    size = ftell(fid);
    if size == -1; error(ferror(fid)); end
    status = fseek(fid, 0, 'bof');
    if status == -1; error(ferror(fid)); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = sf_read_3d_read(fid, pos, size, precision)
% SF_READ_3D_READ   	return data at specified location

status = fseek(fid, pos, 'bof');
[data, cnt] = fread(fid, size, precision);
if cnt ~= prod(size)
    error(['Only read ', num2str(cnt), ' values.']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = sf_read_3d_truncate(data)
% SF_READ_3D_TRUNCATE	truncate data to 8 bits

data(find(data<  0)) =   0;
data(find(data>255)) = 255;
data = uint8(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

