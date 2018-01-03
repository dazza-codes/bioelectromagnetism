function [ avw, machine ] = freesurfer_cor2avw(path,machine)

% freesurfer_cor2avw - Read Freesurfer MRI data (COR-001 to COR-256)
%
% [ avw, machine ] = freesurfer_cor2avw(path, machine)
%
% path - the full path to the COR-??? image files.  If empty,
%        this function uses uigetfile to locate COR-001.
%
% machine - a string, see machineformat in fread for details.
%           The default here is 'ieee-le' but the routine
%           will automatically switch between little and big
%           endian to read any such Analyze header.  It
%           reports the appropriate machine format and can
%           return the machine value.
%
% This function returns an Analyze data structure, with fields:
%
% avw.hdr - a struct with image data parameters.
% avw.img - a 3D matrix of image data (double precision).
%
% See also: avw2cor, avw_read, avw_img_read, avw_hdr_read
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no express or implied warranties
% History:  06/2002, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Create the file header, 256^3, 1 mm^3, uchar, etc
avw = avw_hdr_make;

avw.hdr.dime.dim    =  int16([4 256 256 256 1 0 0 0]);
avw.hdr.dime.pixdim = single([0   1   1   1 0 0 0 0]);


version = '[$Revision: 1.1 $]';
fprintf('\nfreesurfer_cor2avw [v%s]\n',version(12:16));

if ~exist('path','var'),
    fprintf('No input path - see help freesurfer_cor2avw\n');
    [file, path] = uigetfile({'*.*'},'Select COR-001 File');
end
if ~exist('machine','var'), machine = 'ieee-le'; end

[path,name,ext] = fileparts(strcat(path,filesep,'COR-001'));
file = fullfile(path,[name ext]);
fid = fopen(file,'r',machine);

if fid < 0,
    fprintf('Cannot open file %s\n',file);
    fclose(fid);
else
    fclose(fid);
    avw = read_image(avw,path,machine);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ avw ] = read_image(avw,path,machine)

tic;

PixelDim = double(avw.hdr.dime.dim(2));
RowDim   = double(avw.hdr.dime.dim(3));
SliceDim = double(avw.hdr.dime.dim(4));

% Allocate memory for the image
avw.img = zeros(PixelDim,SliceDim,RowDim);

fprintf('...reading ');

% Read one file/slice at a time
for y = 1:SliceDim,

    if y > 1,
        backspaces = repmat('\b',1,22);
    else
        backspaces = '';
    end
    file = sprintf('COR-%03d',y);
    fprintf([backspaces,'%s %s image.'],machine,file);

    [path,name,ext] = fileparts(strcat(path,filesep,file));
    file = fullfile(path,[name ext]);

    % read the whole image into matlab (faster)
    fid = fopen(file,'r',machine); fseek(fid,0,'bof');
    tmp = fread(fid,inf,'uchar=>double');
    fclose(fid);

    % Arrange image into avw.img xyz matrix
    % For Freesurfer COR files (coronal) the voxels are stored with
    % Pixels in 'x' axis (varies fastest) - from patient right to left
    % Rows in   'z' axis                  - from patient superior to inferior
    % Slices in 'y' axis                  - from patient posterior to anterior

    n = 1;
    for z = RowDim:-1:1,
        x = 1:PixelDim;
        avw.img(x,y,z) = tmp(n:n+(PixelDim-1));
        n = n + PixelDim;
    end
end

avw.hdr.hist.orient = uint8(0);

t=toc; fprintf('\n...done (%5.2f sec).\n',t);

return
