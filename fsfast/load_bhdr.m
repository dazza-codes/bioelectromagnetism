function [T, VolRes, VolDim] = load_bhdr(bhdrfile);
%
% [T VolRes VolDim] = load_bhdr(bhdrfile);
%
% Parses the bhdr header file for bshorts/bfloats
%
% T is the affine transform that converts col, row, slice to
% x, y, z in scanner coordinates, assuming that the count starts
% at 0, ie,
%
%      x         col
%      y  = T *  row
%      z        slice
%      1          1
%
% VolRes is a 3x1 vector with the col resolution, row
% resolution, and slice resolution
%
% VolDim is a 3x1 vector with the number of columns, rows, and slices.

if(nargin ~= 1)
  msg = '[T VolRes VolDim] = load_bhdr(bhdrfile)';
  qoe(msg); error(msg);
end

fid = fopen(bhdrfile,'r');
if(fid == -1)
  msg = sprintf('Could not open %s',bhdrfile);
  qoe(msg); error(msg);
end

ok = 1;
while(ok)
  line = fgetl(fid);
  if(line == -1) 
    fclose(fid);
    msg = sprintf('%s appears to be incorrectly formatted\n',bhdrfile);
    qoe(msg); error(msg);
  end
  i = findstr('cols',line) ;
  if(~isempty(i)) break; end
end

                   ncols = sscanf(line,'%*s %d',1);
line = fgetl(fid); nrows = sscanf(line,'%*s %d',1);
line = fgetl(fid); nslices = sscanf(line,'%*s %d',1);

line = fgetl(fid); 
tag = sscanf(line,'%s',1);
if(strmatch(tag,'n_time_points:'))
  ntp = sscanf(line,'%*s %d',1);
  line = fgetl(fid); 
  tag = sscanf(line,'%s',1);
end

slice_res = sscanf(line,'%*s %f',1);

% Get info about the Top Left %
line = fgetl(fid); TLr = sscanf(line,'%*s %f',1);
line = fgetl(fid); TLa = sscanf(line,'%*s %f',1);
line = fgetl(fid); TLs = sscanf(line,'%*s %f',1);
TL = [TLr TLa TLs]'; %'

% Get info about the Top Right %
line = fgetl(fid); TRr = sscanf(line,'%*s %f',1);
line = fgetl(fid); TRa = sscanf(line,'%*s %f',1);
line = fgetl(fid); TRs = sscanf(line,'%*s %f',1);
TR = [TRr TRa TRs]'; %'

% Get info about the Bottom Right %
line = fgetl(fid); BRr = sscanf(line,'%*s %f',1);
line = fgetl(fid); BRa = sscanf(line,'%*s %f',1);
line = fgetl(fid); BRs = sscanf(line,'%*s %f',1);
BR = [BRr BRa BRs]'; %'

% Get the Slice Normal %
line = fgetl(fid); SNr = sscanf(line,'%*s %f',1);
line = fgetl(fid); SNa = sscanf(line,'%*s %f',1);
line = fgetl(fid); SNs = sscanf(line,'%*s %f',1);
Vs = [SNr SNa SNs]'; %'

fclose(fid);

% Compute the Column and Row FOV %
col_fov = sqrt(sum((TR-TL).^2));
row_fov = sqrt(sum((TR-BR).^2));

% Compute the Column and Row Resolution %
col_res = col_fov/ncols;
row_res = row_fov/nrows;

% Compute the Column and Row Normal Vectors %
Vc = (TR-TL)/col_fov;
Vr = (BR-TR)/row_fov;

% The center of the first voxel is TL
VolDim = [ncols nrows nslices]'; %'
VolRes = [col_res row_res slice_res]';%'
P0 = TL;
Mdc = [Vc Vr Vs];
D = diag(VolRes);
T = [Mdc*D P0; 0 0 0 1];

return;
