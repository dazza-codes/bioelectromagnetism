function [LH,RH] = freesurfer_surf_separate(FV, NfacesLH, NvertLH, NfacesRH, NvertRH)

% freesurfer_surf_separate - separate left and right hemi-surfaces
%
% [LH,RH] = freesurfer_surf_separate(FV, NfacesLH, NvertLH, NfacesRH, NvertRH)
% 
% It is assumed that FV input is a struct:
% FV.vertices = Mx3, comprising [LH;RH] vertices (in that order)
% FV.faces = Mx3, comprising [LH;RH] faces (in that order)
% 
% The FV and NfacesLH inputs are essential, all other inputs can be derived
% from these, assuming the LH precedes the RH in the FV struct. This
% assumes that the RH surface was concatenated into FV *after* the LH.
%
% After separating surfaces, try:
% patch('vertices',LH.vertices,'faces',LH.faces,..
%       'facecolor',[1 0 0],'edgecolor','none'); light
% 
% See also freesurfer_read_surf, freesurfer_surf_combine
%


% $Revision: 1.2 $ $Date: 2005/08/25 23:14:06 $

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

% History:  05/2005, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2005/08/25 23:14:06 $';
fprintf('FREESURFER_SURF_SEPARATE [v %s]\n\n',ver(11:15));

if ~exist('FV','var'), 
    error('FV undefined');
end
if isempty(FV),
    error('FV undefined');
end

% check the faces inputs

if ~exist('NfacesLH','var'), 
    error('NfacesLH undefined');
end
if isempty(NfacesLH),
    error('NfacesLH undefined');
end

if ~exist('NfacesRH','var'),
    NfacesRH = length(FV.faces) - NfacesLH;
end
if isempty(NfacesRH),
    NfacesRH = length(FV.faces) - NfacesLH;
end

LHfaceIndex = 1:NfacesLH;
LH.faces = FV.faces(LHfaceIndex,:);
LHvertIndex = unique(LH.faces);
LH.vertices = FV.vertices(LHvertIndex,:);

% check the vertex inputs

if ~exist('NvertLH','var'),
    NvertLH = length(LHvertIndex);
end
if isempty(NvertLH),
    NvertLH = length(LHvertIndex);
end

if ~exist('NvertRH','var'),
    NvertRH = length(FV.vertices) - NvertLH;
end
if isempty(NvertRH),
    NvertRH = length(FV.vertices) - NvertLH;
end

if ~isequal(LHvertIndex, [1:NvertLH]'),
    warning('LH faces do not index [1:NvertLH]');
end
if (NvertLH+NvertRH) ~= length(FV.vertices),
    error('(NvertLH+NvertRH) ~= length(FV.vertices)');
end
if (NfacesLH+NfacesRH) ~= length(FV.faces),
    error('(NfacesLH+NfacesRH) ~= length(FV.faces)');
end

RHvertIndex = (NvertLH+1):(NvertLH+NvertRH);
RH.vertices = FV.vertices(RHvertIndex,:);
RHfaceIndex = (NfacesLH+1):(NfacesLH+NfacesRH);
RH.faces = FV.faces(RHfaceIndex,:) - length(LH.vertices);

return
