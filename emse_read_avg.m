function [EMSE] = emse_read_avg(filename)

% emse_read_avg - Load EMSE .avg data (actually ascii format)
% 
% Useage: [EMSE] = emse_read_avg(filename)
%
% where 'filename' is the full path + fileprefix + filextension
%
% The returned struct has the following fields:
% 
% EMSE.channels
% EMSE.pnts
% EMSE.rate       - sample rate (msec)
% EMSE.xmin       - prestim baseline period (msec)
% EMSE.volt       - potential floating point matrix, 
%                   size [points,channels]
% 
% No variance data is yet read or returned
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2000, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('EMSE_READ_AVG [v %s]\n',ver(11:15));

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

if exist(file) ~= 2,
  lookfile = which(file);
  if isempty(lookfile),
    msg = sprintf('Cannot locate %s\n', filename);
    error(msg);
  else
    file = lookfile;
  end
end

fprintf('...reading: %s\n', file);

fid = fopen(file);

version   = fscanf(fid,'%d',1);
file_type = fscanf(fid,'%d',1);
minor_rev = fscanf(fid,'%d',1);

if isempty(version),
  EMSE.channels = [];
  EMSE.pnts = [];
  EMSE.rate = [];
  EMSE.xmin = [];
  EMSE.volt = [];
  fprintf('...this is not an EMSE file.\n...it might be a Neuroscan file.\n');
  return
end

fprintf('...Version = %d, File-Type = %d, Minor_Revision = %d\n',...
  version,file_type,minor_rev);

unknown  = fscanf(fid,'%d',1);
channels = fscanf(fid,'%d',1);
points   = fscanf(fid,'%d',1);
samples  = fscanf(fid,'%f',1) * 1000; % msec sample rate
unknown  = fscanf(fid,'%f',1);
baseline = fscanf(fid,'%f',1) * -1000; % msec baseline
unknown  = fscanf(fid,'%d',1);
unknown  = fscanf(fid,'%d',1);


fprintf('...Sample Rate (msec) = %6.3f, Baseline (msec) = %6.3f\n',...
  samples,baseline);

for i = 1:channels,
  discard = fscanf(fid,'%d',2)';
end

volt = zeros(points,channels);

for i = 1:points,
  volt(i,:) = fscanf(fid,'%f',channels)';
end

fclose(fid);

fprintf('...Points (rows) = %d, Channels (cols) = %d\n',points,channels);

EMSE.channels = channels;
EMSE.pnts = points;
EMSE.rate = samples;
EMSE.xmin = baseline;
EMSE.volt = volt;

return
