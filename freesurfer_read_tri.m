function [vertices,faces] = freesurfer_read_tri(file)

% freesurfer_read_tri - Read FreeSurfer tesselation (.tri)
% 
% [vertices,faces] = freesurfer_read_tri(file)
% 
% This function will load an ascii file that contains a one
% line specification of the number of vertices followed 
% by rows of vertex points.  It then reads a one line
% specification of the number of faces followed by rows
% of face indices into the vertex rows.  Each vertex row 
% contains a vertex index number and 3 x,y,z coordinates. 
% Each face row contains a face index and three vertex 
% indices.  Vertices in the .tri file are indexed from one
% and those returned are indexed from one.
% 
% See also the mesh_freesurfer2matlab function to load
% the tesselations that are created by the mris_convert 
% function of freesurfer, which have a different text format 
% from those of the BEM .tri files.
% 
% The freesurfer tesselations may contain too many faces
% for efficient computations.  If so, try 'reducepatch'.
%
% The vertex coordinates are in mm.  The FreeSurfer coordinate
% system is too confusing to explain here; see the FreeSurfer
% homepage or some nice documentation by Graham Wideman.
%
% The returned matrices can be input to the patch command, like so:
%
% Hpatch = patch('Vertices',vertices,'Faces',faces,...
%                'EdgeColor',[.8 .8 .8],'FaceColor',[0.9 0.9 0.9]);
%
% This will plot the mesh as a patch object.  See the patch command
% and matlab help for more information on coloring this object.
%


% Copyright (C) 2002  Darren L. Weber
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

% History:  03/2002 Darren.Weber_at_radiology.ucsf.edu
%           02/2004 Darren.Weber_at_radiology.ucsf.edu
%                   previously called mesh_freesurferTRI2matlab
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.2 $ $Date: 2005/01/21 04:33:50 $';
fprintf('FREESURFER_READ_TRI [v %s]\n',ver(11:15));

fid = fopen(file,'r');

if isequal(fid,-1),
    msg = sprintf('Could not open file: %s',file);
    error(msg);
else
    
    fprintf('...reading FreeSurfer triangulation (.tri)\n');
    
    tic;
    
    % Check for comment on first line of file
    frewind(fid); temp = fscanf(fid,'%s',1); frewind(fid);
    if findstr(temp,'#'), temp = fgetl(fid); end
    
    % Read vertices
    Nvertices = fscanf(fid,'%d',1);
    fprintf('...Reading %d Vertices\n',Nvertices);
    vertices = fscanf(fid,'%f',[4,Nvertices]);
    % remove first row (index) and translate
    vertices = vertices(2:4,:)';
    
    % Read faces
    Nfaces = fscanf(fid,'%d',1);
    fprintf('...Reading %d Faces\n',Nfaces);
    faces = fscanf(fid,'%d',[4,Nfaces]);
    % remove first row (index) & translate
    faces = faces(2:4,:)';
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n\n',t);
    
end

return
