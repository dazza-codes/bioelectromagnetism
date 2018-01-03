function [p] = cnt_open(p,parent)

% cnt_open - function to handle various eeg_load commands
% 
% Usage: [p] = cnt_open(p,[parentgui])
% 
% p is a parameter structure. See eeg_toolbox_defaults for more 
% information on this parameter structure.
% 
% In this function, p must contain the fields:
% 
% p.cnt.path - the directory location of the file to load
% p.cnt.file - the name of the file to load
% p.cnt.type - the file format string, one of:
% 
% 'Scan4x'
% 'Scan3x'
%
% These are the only CNT file types currently supported. 
% See functions eeg_load* for details.
% 
% The most important return value is the CNT data in 
% p.cnt.  If the file format is scan4x, various 
% parameters are returned also.
% 
% See also: EEG_LOAD, EEG_LOAD_ASCII, EMSE_READ_AVG, 
%           EEG_LOAD_SCAN4_CNT, EEG_LOAD_SCAN3_CNT
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no express or implied warranties
% History:  10/2003, Darren.Weber_at_radiology.ucsf.edu
%                    created
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('p','var'),[p] = eeg_toolbox_defaults; end

eegversion = '$Revision: 1.1 $';
fprintf('CNT_OPEN [v %s]\n',eegversion(11:15)); tic;


[path,name,ext] = fileparts(strcat(p.cnt.path, filesep, p.cnt.file));
file = fullfile(path,[name ext]);

if ~isequal(exist(file),2),
  lookfile = which(file);
  if isempty(lookfile),
    msg = sprintf('...cannot locate %s\n', file);
    error(msg);
  else
    file = lookfile;
  end
end

type = lower(p.cnt.type)

switch type,
  
  case 'scan4x_cnt',
    
    p.cnt = eeg_load_scan4_cnt(file);
    
  case 'scan3x_cnt',
    
      msg = sprintf('...failed to load scan3.x datafile: %s\n',file);
      error(msg);
      
  otherwise,
    msg = sprintf('\nPlease specify data type: Scan4x_cnt | Scan3x_cnt\n\n');
    error(msg);
end


t = toc; fprintf('...done (%6.2f sec).\n\n',t);

return
