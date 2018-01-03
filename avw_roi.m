function [roi] = avw_roi(avw,position,shape)

% avw_roi - extract a region of interest from avw.img
%
% [roi] = avw_roi(avw,[position],[shape])
%
% avw - the analyze struct with avw.hdr and avw.img
%
% position - the centroid voxel for the roi, the default is the volume
% center.  For example, 256^3 voxels would have a center at [128,128,128].
%
% shape - the shape of the roi.  This is a struct:
%
%            shape.type = 'block'
%            shape.size = [X,Y,Z] voxels
%
% The default is 'block' and the whole volume.  Note that the block type
% encompasses any 3D square or rectangle.  Further development may include
% sphere and ellipse shapes, with mm parameters.
%
% roi.index - The volume indices into avw.img that comprise the region of
% interest (roi).
% roi.value - The values at the voxels in avw.img that comprise the roi.
% The returned value can be used with avw_stats.  The default is the entire
% volume, ie: roi.value = avw.img;
%
% see also:  avw_stats
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



version = '[$Revision: 1.1 $]';
fprintf('\nAVW_ROI [v%s]\n',version(12:16));  tic;

if ~exist('avw','var'),
  doc avw_roi;
  msg = sprintf('...no input avw\n');
  error(msg);
elseif isempty(avw),
  doc avw_roi;
  msg = sprintf('...empty input avw\n');
  error(msg);
end

if ~exist('shape','var'),
  fprintf('...no input shape, using the whole volume\n');
  roi.value = avw.img;
  return
else
  if ~isfield(shape,'type'),
    fprintf('...no input shape.type, using block\n');
    shape.type = 'block';
  end
  if ~isfield(shape,'size'),
    fprintf('...no input shape.size, using [5,5,5]\n');
    shape.size = [ 5,5,5 ];
  end
end

s = size(avw.img);
xdim = s(1);
ydim = s(2);
if length(s) > 2, zdim = s(3);
else              zdim = 1;
end

if ~exist('position','var'),
  fprintf('...no input position, using volume center\n');
  position = [1,1,1];
  if xdim > 1, position(1) = floor(xdim/2); end
  if ydim > 1, position(2) = floor(ydim/2); end
  if zdim > 1, position(3) = floor(zdim/2); end
end



switch shape.type,
  
  case 'block',
    
    fprintf('...defining [%d,%d,%d] ''block'' region of interest.\n',...
      shape.size(1),shape.size(2),shape.size(3));
    
    xrange = floor(shape.size(1) / 2);
    yrange = floor(shape.size(2) / 2);
    zrange = floor(shape.size(3) / 2);
    
    xroi = [ position(1) - xrange, position(1) + xrange ];
    if xroi(1) < 1, xroi(1) = 1; end
    if xroi(2) > xdim, xroi(2) = xdim; end
    roi_xindex = xroi(1):xroi(2);
    
    yroi = [ position(2) - yrange, position(2) + yrange ];
    if yroi(1) < 1, yroi(1) = 1; end
    if yroi(2) > ydim, yroi(2) = ydim; end
    roi_yindex = yroi(1):yroi(2);
    
    zroi = [ position(3) - zrange, position(3) + zrange ];
    if zroi(1) < 1, zroi(1) = 1; end
    if zroi(2) > zdim, zroi(2) = zdim; end
    roi_zindex = zroi(1):zroi(2);
    
    roi.index = [roi_xindex;roi_yindex;roi_zindex];
    roi.value = avw.img(roi_xindex,roi_yindex,roi_zindex);
    
  case 'sphere',
    
    %fprintf('...defining [%d,%d,%d] ''sphere'' region of interest.\n',...
    %  position(1),position(2),position(3));
    
    fprintf('...sorry, sphere shape calculations are not yet available\n');
    roi.index = [];
    roi.value = [];
    
    %Rroi = shape.size(1);
    
    
  otherwise,
    
    msg = sprintf('...no support for shape.type = %s',shape.type);
    error(msg);
    
end

return
