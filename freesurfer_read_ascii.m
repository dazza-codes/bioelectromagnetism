function [vertices,faces] = freesurfer_read_ascii(file)

% freesurfer_read_ascii - Read FreeSurfer tesselation (.txt)
%
% [vertices,faces] = freesurfer_read_ascii(file)
%
% This function will load an ascii file that contains a one
% line specification of the number of vertices and faces,
% followed by rows of vertex points and then rows of face
% indices into the vertex rows.  Each vertex row contains
% 3 x,y,z coordinates and a colour/potential value.  Each 
% face row contains three vertex indices followed by a
% colour/potential value.  Vertices in the .txt file are
% indexed from zero, but those returned are indexed from
% one.
% 
% See also the freesurfer_read_tri function to load
% the data from the BEM .tri tesselations, which have a
% different text format from those of the .txt/.asc files 
% that are created by the mris_convert function of freesurfer.
% 
% The freesurfer tesselations may contain too many faces
% for efficient computations.  If so, try 'reducepatch'.
% 
% The vertex coordinates are in mm.  The FreeSurfer coordinate
% system is too confusing to explain here; see the FreeSurfer
% homepage or some nice documentation by Graham Wideman.
%
% The returned values can be input to the patch command, like so:
%
%    Hpatch = patch('Vertices',vertices,'Faces',faces,...
%                   'EdgeColor',[.8 .8 .8],'FaceColor',[0.9 0.9 0.9]);
%
% This will plot the mesh as a patch object.  See the patch command
% and matlab help for more information on coloring this object.
%
% Freesurfer website: http://surfer.nmr.mgh.harvard.edu/
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:34 $

% Licence:  GNU GPL, no implied or express warranties
% History:  03/2002 Darren.Weber_at_radiology.ucsf.edu
%           02/2004 Darren.Weber_at_radiology.ucsf.edu
%                   previously called mesh_freesurfer2matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(file,'r');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    fprintf('...Reading FreeSurfer Tesselation Data\n');
    
    tic;
    
    % Check for comment on first line of file
    frewind(fid); temp = fscanf(fid,'%s',1); frewind(fid);
    if findstr(temp,'#'), temp = fgetl(fid); end
    
    % Read number of vertices/faces
    Nvertices = fscanf(fid,'%d',1);
    Nfaces    = fscanf(fid,'%d',1);
    
    % Read vertices
    fprintf('...Reading %d Vertices\n',Nvertices);
    vertices = fscanf(fid,'%f',[4,Nvertices]);
    % remove last row (all zeros) and translate
    vertices = vertices(1:3,:)';
    
    % Read faces
    fprintf('...Reading %d Faces\n',Nfaces);
    faces = fscanf(fid,'%d',[4,Nfaces]);
    % remove last row (all zeros), translate and
    % add 1 because FreeSurfer vertices start at zero
    faces = faces(1:3,:)' + 1;
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n',t);
    
end

return
