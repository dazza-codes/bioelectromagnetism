function fname = TalDaemon_run(fname)

% fname = TalDaemon_run(fname)
%
% The function assumes the Talairach Daemon command line client is
% installed into a directory called 'TalairachDaemon' below the path where
% this function is installed.  It runs the daemon on the input file and
% returns the filename of the output (*.td) file.
%
% See online: http://ric.uthscsa.edu/TDinfo
%             http://www.brainmap.org/
%

% $Revision: 1.2 $ $Date: 2007/11/26 23:04:14 $

% Copyright (C) 2005  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% History:  08/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2007/11/26 23:04:14 $';
fprintf('TALDAEMON_RUN [v %s]\n',ver(11:15));

if(nargin < 1)
  help TalDaemon_run;
  return;
end

TalDaemonPath = TalDaemon_tools_path;

% The following is not sufficient to identify the right client programs,
% because some of them are compiled for different flavours of unix.
if ispc,
  TalDaemonExe = fullfile(TalDaemonPath,'excel2td.exe');
else
  TalDaemonExe = fullfile(TalDaemonPath,'tdc');
end

if ~exist(TalDaemonExe),
  msg = sprintf('cannot locate Talairach Daemon Client,\n%s\n',TalDaemonExe);
end

cwd = pwd; % store the pwd to come back to it
[TDfilepath,TDfile,TDext] = fileparts(fname);
cd(TDfilepath)
if ispc,
  commandstr = sprintf('%s 2, %s', TalDaemonExe, [TDfile,TDext]);
else
  commandstr = sprintf('./%s 2, %s', TalDaemonExe, [TDfile,TDext]);
end
[s,w] = system(commandstr);

if s < 0,
  error(w);
end

cd(cwd);

fname = fullfile(TDfilepath,[TDfile,TDext,'.td']);

return
