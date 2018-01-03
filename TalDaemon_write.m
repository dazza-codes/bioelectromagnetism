function TalDaemon_write(fname,TALdata)

% TalDaemon_write(fname,TALdata)
%
% fname is a full path to a text file
% TALdata is integer or float values, Nx3 XYZ Talairach coordinates; for
% processing by the Talairach Daemon, these values are rounded off to
% integers and output in N lines of tab delimited ascii text.
%
% See: TalDaemon_run, TalDaemon_read
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $ $Date: 2005/09/10 05:00:54 $';
fprintf('TALDAEMON_WRITE [v %s]\n',ver(11:15));

fid = fopen(fname, 'w');
if(fid == -1)
    msg = sprintf('Could not open %s\n',fname);
    error(msg)
end

% convert TALdata to integers
TALint = round(TALdata);

[s1,s2] = size(TALint);

if s2 == 3,
    count = fprintf(fid,'%d\t%d\t%d\n',TALint');
else
    count = fprintf(fid,'%d\t%d\t%d\n',TALint);
end

fclose(fid);

return
