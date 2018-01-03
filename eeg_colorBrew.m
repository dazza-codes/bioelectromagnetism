function [p] = eeg_colorBrew(p,cname,ncolors)

% eeg_colorBrew - import colormaps from colorBrew
%
%   [p] = eeg_colorBrew(p,cname,ncolors)
%
%   This function calls colorBrew.  See colorBrew for explanation of cname
%   and ncolors.  Explore colorBrew first, if unsure about what to specify
%   for this function.
%
%   Returns a colormap matrix [ncolors, 3] into p.colorMap.map, which can be
%   used with the matlab colormap command.  For example:
%
%   [p] = eeg_colorBrew(p,'Spectral',11);
%   colormap(p.colorMap.map);
%
%   The default p struct is created by eeg_toolbox_defaults.  It contains
%   p.colorMap and this function imports a colorBrew map into the
%   p.colorMap.map array.  The defaults are:
%   
%   style:  'Spectral'   (cname input to colorBrew)
%   length: 11           (ncolors input to colorBrew)
%   exp:    0            (has no meaning or purpose in this function)
%   Cmin:   0            (always 0, colorBrew returns values from 0-1)
%   Cmax:   1            (always 1, colorBrew returns values from 0-1)
%   plot:   0            (no plot)
%

% $Revision: 1.2 $ $Date: 2007/04/30 18:00:49 $

% Licence:  GNU GPL, no implied or express warranties
% History:  02/2007, Darren.Weber_at_radiology.ucsf.edu, created
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),
  error('EEG_COLORBREW: No input p struct, use eeg_toolbox_defaults.');
elseif strmatch(class(p), 'struct', 'exact'),
  if ~isfield(p,'colorMap'), p.colorMap = []; end
  if ~isfield(p.colorMap,'style'),  p.colorMap.style   = 'Spectral'; end
  if ~isfield(p.colorMap,'length'), p.colorMap.length  = 11;  end
  if ~isfield(p.colorMap,'exp'),    p.colorMap.exp     = 0;   end
  if ~isfield(p.colorMap,'Cmin'),   p.colorMap.Cmin    = 0;   end
  if ~isfield(p.colorMap,'Cmax'),   p.colorMap.Cmax    = 1;   end
  if ~isfield(p.colorMap,'plot'),   p.colorMap.plot    = 0;   end
else
  msg = sprintf('p is unknown type: %s', class(p));
  error(msg);
end

% Ensure that Length is defined
if ~isfield(p.colorMap,'length'),
  p.colorMap.length  = 6;
else
  p.colorMap.length = fix(p.colorMap.length);
end

if ~exist('cname','var'),
  cname = 'Spectral';
end
if isempty(cname),
  cname = 'Spectral';
end

if ~exist('ncolors','var'),
  ncolors = 11;
end
if isempty(ncolors),
  ncolors = 11;
end

p.colorMap.style = cname;
p.colorMap.length = ncolors;
p.colorMap.exp     = 0;
p.colorMap.Cmin    = 0;
p.colorMap.Cmax    = 1;
p.colorMap.plot    = 0;

%A color map matrix may have any number of rows, but it must have
%exactly 3 columns.  Each row is interpreted as a color, with the
%first element specifying the intensity of red light, the second 
%green, and the third blue.  Color intensity can be specified on the
%interval 0.0 to 1.0.

cb = colorBrew(cname,ncolors);

p.colorMap.map = flipud(cb.colorRGB);


% Plot colormap
if (p.colorMap.plot),
  figure('NumberTitle','off','Name','EEG Toolbox Colormap');
  colormap(p.colorMap.map);
  rgbplot(p.colorMap.map);
  colorbar
  axis tight
end

return
