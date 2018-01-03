function TDdata = TalDaemon_read(fname)

% TDdata = TalDaemon_read(fname)
%
% returns a structure of the fields from the Talairach daemon report
% TDdata.xyz
% TDdata.hemi
% TDdata.lobe
% TDdata.gyrus
% TDdata.matter
% TDdata.Brodmann
%
% However, the Talairach Daemon does not have constant output formats.  The
% command line version and the GUI client have different output formats (on
% windows at least).  Hence, this function attempts to detect and adapt to
% the different formats, but it is not a fail proof routine!
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
fprintf('TALDAEMON_READ [v %s]\n',ver(11:15));

fid = fopen(fname,'r');
if(fid == -1)
    msg = sprintf('Could not open %s\n',fname);
    error(msg)
end

fprintf('...reading %s\n',fname);

% If the Talairach Daemon had constant output formats, the following would
% not be necessary.  However, the command line version and the GUI client
% have different output formats (on windows at least).  Hence, the
% following mess is an attempt to work with this.  It is not a fail proof
% detection and parsing routine!

str = fgetl(fid);
if str == -1,
    % file is empty
    error('file is empty');
end
if strfind(str,','),
    % this file is comma separated values (csv)
    delimit = ',';
    
    %eg,
    % str =
    % 0, -7, 22, 35, Left Cerebrum,Frontal Lobe,Cingulate Gyrus,Gray Matter,Brodmann area 32,Range=1
    % Field1 = index
    % Field2-4 = xyz
    
    % there may be a way to automatically determine the length of this
    % format string by parsing 'str'.  By manual inspection, this format
    % contains 4 leading integers, starting with an entry index of zero.
    % It then contains a series of comma separated strings
    fieldLoc = strfind(str,delimit);
    fieldN = length(fieldLoc) + 1;
    Nindex = 1:4;
    format = '';
    for i = 1:Nindex(end), format = [format,'%d']; end
    for i = 4:fieldN, format = [format,'%s']; end
elseif strfind(str,'\t'),
    % this file is tab delimited
    delimit = '\t';
    fieldLoc = strfind(str,delimit);
    fieldN = length(fieldLoc) + 1;
    Nindex = 1:3;
    format = '';
    for i = 1:Nindex(end), format = [format,'%d']; end
    for i = 3:fieldN, format = [format,'%s']; end
else
    error('unknown file format, not comma or tab delimited');
end


fseek(fid,0,-1);
switch delimit,
    case ',',
        C = textscan(fid,format,'delimiter',delimit);
    case '\t',
        C = textscan(fid,format,'delimiter',delimit);
end
fclose(fid);

N = length(Nindex);
if N > 3,
    TDdata.xyz = [ C{2} C{3} C{4} ];
else
    TDdata.xyz = [ C{1} C{2} C{3} ];
end
TDdata.hemi = C{N+1};
TDdata.lobe = C{N+2};
TDdata.gyrus = C{N+3};
TDdata.matter = C{N+4};

% convert BA values to numbers
BA = strrep(C{N+5},'Brodmann area ','');
BA = str2double(BA);
TDdata.Brodmann = BA;

return
