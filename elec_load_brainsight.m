function [brainsight] = elec_load_brainsight(filename)

% elec_load_brainsight - read ascii export of BrainSight electrode data
% 
% Load an ascii electrode file and return the electrode labels and
% coordinates.  This function can read and return Cartesian (x,y,z) or
% spherical (theta,phi,r) coordinates.
%
% [brainsight] = elec_load_brainsight(filename)
% 
% where:
% 
% file = '<path><filename>'
% 
% The file format is that of BrainSight ascii export files.  The file
% contains a short header and then two data sections, first a list of x,y,z
% locations for each electrode and then a second listing, according to
% electrode label.  Each data field is separated by spaces.  For example,
% 
% fz
% X.X Y.Y Z.Z
% 
% Example result:
%
% brainsight = 
%
%        n: 72
%    units: 'mm'
%    label: {72x1 cell}
%        x: [72x1 double]
%        y: [72x1 double]
%        z: [72x1 double]
%
% All coordinates are in brainsight.units (mm).
% 

% $Revision: 1.1 $ $Date: 2007/02/14 07:17:32 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2007, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $';
fprintf('\nELEC_LOAD_BRAINSIGHT [v %s]\n',ver(11:15));

tic;

[path,name,ext] = fileparts(filename);
file = fullfile(path,[name ext]);

fprintf('...loading electrodes from:\n\t%s\n', file);

brainsight = read_brainsight(file);

fprintf('...loaded %d electrodes\n', size(brainsight.x,1));

t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function brainsight = read_brainsight(file),

fid = fopen(file);
if fid < 0,
  msg = sprintf('cannot open file: %s\n',file);
  error(msg);
end


% ------------------------------------------------
% Read the short header, eg:
%%# Electrodes
%%NumberPositions=  72
%%UnitPosition mm

tmp = fgetl(fid);
if strmatch('# Electrodes', tmp, 'exact') == 1,

  % read the number of positions or fail
  tmp = fgetl(fid);
  if strmatch('NumberPositions=', tmp) == 1,
    brainsight.n = str2num(strrep(tmp,'NumberPositions=',''));
  else
    error('failed to find ''NumberPositions='' in header');
  end
  
  % read the position units or fail
  tmp = fgetl(fid);
  if strmatch('UnitPosition', tmp) == 1,
    brainsight.units = strrep(tmp,'UnitPosition ','');
  else
    error('failed to find ''UnitPosition'' in header');
  end
  
else
  error('failed to read header correctly');
end


% ------------------------------------------------
% initialize data fields

brainsight.label = cell(brainsight.n, 1);
brainsight.x = zeros(brainsight.n, 1);
brainsight.y = brainsight.x;
brainsight.z = brainsight.x;
%brainsight.hsp = [];  % for head-shape points


% ------------------------------------------------
% read the position x,y,z values

fprintf('...reading positions x y z values\n');

tmp = fgetl(fid);
if strmatch('Positions', tmp) ~= 1,
  error('failed to find ''Positions'' in header');
end

for posN = 1:brainsight.n,
  
  % read the next line of the file
  data = fgetl(fid);
  if data < 0,                           % we are at EOF
    break
  end
  
  data = str2num(data);
  
  brainsight.x(posN) = data(1);
  brainsight.y(posN) = data(2);
  brainsight.z(posN) = data(3);
  
  %fprintf('line %4d: %g %g %g\n', posN, data)
  
end

if posN < brainsight.n,
  msg = sprintf('did not read %d x,y,z values for positions', brainsight.n);
  warning(msg)
end

fprintf('...read %d positions\n', posN)



% -------------------------------------------------------------------------
% for some strange reason, their data repeats all over again, this time
% using electrode labels, so now we read the position label and confirm the
% x,y,z values (we could just skip the entire read above and just do it once
% here, but let's be anal retentive about this relatively new data format
% for now).

for posN = 1:brainsight.n,
  
  % read the next line of the file
  label = fgetl(fid);
  if label < 0,                           % we are at EOF
    break
  end
  
  % the first line is the electrode label (lower case)
  brainsight.label(posN) = {label};
  
  
  % read the next line of the file
  data = fgetl(fid);
  if data < 0,                           % we are at EOF
    break
  end
  
  % run some anal retentive data verification, but this might be susceptible
  % to rounding errors for floats
  data = str2num(data);
  if brainsight.x(posN) ~= data(1),
    msg = sprintf('inconsistent data specification, at position %d x', posN);
    warning(msg);
  end
  if brainsight.y(posN) ~= data(2),
    msg = sprintf('inconsistent data specification, at position %d y', posN);
    warning(msg);
  end
  if brainsight.z(posN) ~= data(3),
    msg = sprintf('inconsistent data specification, at position %d z', posN);
    warning(msg);
  end
  
  % Oh - the data is inconsistent!  The error is within the first decimal
  % place, so this is something more than rounding errors.  (DLW 02/13/2007)
  
  
  % We should not do the following without knowing why, but for the sake of
  % convenience, let's use the last set of values for now.  There is no
  % assumption here about what set of values are 'correct'.
  brainsight.x(posN) = data(1);
  brainsight.y(posN) = data(2);
  brainsight.z(posN) = data(3);
  
  
end


fclose(fid);

return




% $$$ fprintf('...searching for fiducials...');
% $$$ 
% $$$ n = 0;
% $$$ while n < 5,
% $$$   tmp = fgetl(fid);
% $$$   if tmp < 0, break; end
% $$$   
% $$$   tmp = lower(tmp);
% $$$   
% $$$   if strfind(tmp,'nasion'),
% $$$     tmp = sscanf(tmp,'%s %d %f %f %f');
% $$$     brainsight.nasion = [tmp(end-2) tmp(end-1) tmp(end)];
% $$$     n = n + 1;
% $$$     continue;
% $$$   end
% $$$   if strfind(tmp,'left'),
% $$$     tmp = sscanf(tmp,'%s %d %f %f %f');
% $$$     brainsight.lpa = [tmp(end-2) tmp(end-1) tmp(end)];
% $$$     n = n + 1;
% $$$     continue;
% $$$   end
% $$$   if strfind(tmp,'right'),
% $$$     tmp = sscanf(tmp,'%s %d %f %f %f');
% $$$     brainsight.rpa = [tmp(end-2) tmp(end-1) tmp(end)];
% $$$     n = n + 1;
% $$$     continue;
% $$$   end
% $$$   if strfind(tmp,'centroid'),
% $$$     tmp = sscanf(tmp,'%s %d %f %f %f');
% $$$     brainsight.origin = [tmp(end-2) tmp(end-1) tmp(end)];
% $$$     n = n + 1;
% $$$     continue;
% $$$   end
% $$$   if strfind(tmp,'ref'),
% $$$     tmp = sscanf(tmp,'%s %d %f %f %f');
% $$$     brainsight.ref = [tmp(end-2) tmp(end-1) tmp(end)];
% $$$     n = n + 1;
% $$$     continue;
% $$$   end
% $$$ end
% $$$ 
% $$$ frewind(fid);
% $$$ fprintf('done\n');
