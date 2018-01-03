function [Channel] = emse_elp2brainstorm(elp,chanFile)

% emse_elp2brainstorm - Convert EMSE elp to brainstorm channel file
% 
% The EMSE elp struct is returned from emse_read_elp.  The elp data
% structure is converted into the brainstorm format and returned.
% 
% Useage: Channel = emse_elp2brainstorm(elp,[brainstormChanFile])
% 
% elp = see emse_read_elp for more details
% 
% brainstormChanFile = a full path to a channel.mat file.  If this is
% empty, the function will not save an output file.
%
% Channel is an array of structures.  The fields are:
% 
%  Loc     - a 3x2 matrix of electrode and reference coordinates.  Each
%            column contains [X;Y;Z] values.
%  Orient  - a corresponding matrix of sensor orientations for MEG; 
%             all zero for EEG.
%  Weight  - a vector of relative or absolute weights (eg, gain);
%            all ones for this routine.
%  Type    - a character string, 'EEG' in this function.
%  Name    - a charater string indicating the electrode name.
%  Comment - a charater string indicating the reference electrode. Empty
%            for active electrodes and 'EEG REF' for the reference.
% 
% See brainstorm website at http://neuroimage.usc.edu/, including a
% download pdf file describing the brainstorm database formats.
% 

% $Revision: 1.1 $ $Date: 2007/05/08 21:38:49 $

% Licence:  GNU GPL, no express or implied warranties
% History:  05/2007, Darren.Weber_at_radiology.ucsf.edu
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('elp', 'var'),
  error('no input elp struct')
end
if isempty(elp),
  error('empty elp struct')
end

if ~exist('chanFile', 'var'),
  chanFile = '';
end

tic;

  ver = '$Revision: 1.1 $';
  fprintf('\nEMSE_ELP2BRAINSTORM [v %s]\n',ver(11:15));
  fprintf('...Converting to brainstorm structure.\n');

  for i=1:length(elp.x),
    Channel(i).Loc = [[elp.x(i) elp.y(i) elp.z(i)]',elp.ref'];
    Channel(i).Orient = [];     % used for MEG rather than EEG
    Channel(i).Weight = 1;      % Like Amplification
    Channel(i).Type = 'EEG';
    Channel(i).Name = elp.name{i};
    Channel(i).Comment = '';
  end
  Channel(i+1).Loc = [elp.ref',elp.ref'];
  Channel(i+1).Orient = [];
  Channel(i+1).Weight = 1;
  Channel(i+1).Type = 'EEG';
  Channel(i+1).Name = 'EEG REF';
  Channel(i+1).Comment = 'EEG REF';

  if ~isempty(chanFile),
    fprintf('...saving BrainStorm channel data to:\n...%s\n',chanFile);
    save(chanFile, 'Channel');
  end

  t = toc; fprintf('...done (%6.2f sec).\n\n',t);

  return
