function meshtype = freesurfer_meshtype(file),

% freesurfer_meshtype - get FreeSurfer mesh type from file name
%
% meshtype = freesurfer_meshtype(file),
%
% meshtype might be any of these string values:
%
% 'orig'
% 'smoothwm'
% 'white'
% 'pial'
% 'brain' 'inner skull' - brain actually approximates the inner skull
% 'scalp'
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no express or implied warranties
% History:  02/2004 Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  if findstr(file,'smoothwm'),  meshtype = 'smoothwm';
  elseif findstr(file,'orig'),  meshtype = 'orig';
  elseif findstr(file,'white'), meshtype = 'white';
  elseif findstr(file,'pial'),  meshtype = 'pial';
  elseif findstr(file,'brain'), meshtype = 'inner skull';
  elseif findstr(file,'scalp'), meshtype = 'scalp';
  elseif findstr(file,'head'),  meshtype = 'scalp';
  else,                         meshtype = 'unknown';
  end
  
return
  
