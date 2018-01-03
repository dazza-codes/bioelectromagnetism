function [TalDaemonPath] = TalDaemon_tools_path

% TalDaemon_tools_path - locate the Talairach Daemon installation path
%
% [TalDaemonPath] = TalDaemon_tools_path
%
% This function looks for the installation path of
% TalDaemon_run.m in the matlab path.  If not 
% found, update the matlab path (try addpath). It
% assumes the Talairach daemon command line client
% is installed into a directory below the path to
% TalDaemon_run.m
%


% $Revision: 1.1 $ $Date: 2005/09/10 05:00:54 $

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TDtoolsPath = fileparts(which('TalDaemon_run'));
if isempty(TDtoolsPath),
    msg = sprintf('Cannot find TalDaemon_run on the matlab path.\nTry the addpath command.\n\n');
    error(msg);
end

% assume the Talairach Daemon client is installed below this path
TalDaemonPath = fullfile(TDtoolsPath,'TalairachDaemon','');
addpath(TalDaemonPath,'-end');

return
