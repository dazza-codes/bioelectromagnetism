function [vertices,faces] = geomview_read_off(file)

% geomview_read_off - Read GeomView .off mesh file format
%
% [vertices,faces] = geomview_read_off(file)
%
% This function will load an ascii file that contains a one line
% specification of the 'OFF' file format followed by another line to
% specify number of vertices, faces and edges, followed by rows of vertex
% points and then rows of face indices into the vertex rows.  Each vertex
% row contains x,y,z coordinates.  Each face row contains the number of
% vertices in the face and then the vertex indices.  Vertices in the .off
% file are indexed from zero, but those returned are indexed from
% one.
% 
% See http://www.geomview.org/docs/oogltour.html
%
% The returned values can be input to the patch command, like so:
%
% Hpatch = patch('Vertices',vertices,'Faces',faces,...
%                'EdgeColor',[.8 .8 .8],'FaceColor',[0.9 0.9 0.9]);
%
% This will plot the mesh as a patch object.  See the patch command
% and matlab help for more information on coloring this object.
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Licence:  GNU GPL, no implied or express warranties
% History:  11/2004 Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(file,'r');

if isequal(fid,-1),
    S=sprintf('Could not open file: "%s"',file);
    error(S);
else
    
    fprintf('...Reading GeomView .off mesh file\n');
    
    tic;
    
    % Check for comment on first line of file
    frewind(fid); temp = fgetl(fid); frewind(fid);
    if findstr(temp,'OFF'),
        % good, we have the right format file
        temp = fgetl(fid);
    else
        error('This is not a GeomView .off file, the first line must contain the header ''OFF''');
    end
    
    % Read number of vertices/faces/edges
    Nvertices = fscanf(fid,'%d',1);
    Nfaces    = fscanf(fid,'%d',1);
    Nedges    = fscanf(fid,'%d',1);
        
    % Read vertices
    fprintf('...Reading %d Vertices\n',Nvertices);
    vertices = fscanf(fid,'%f',[3,Nvertices]);
    % select first 3 rows and translate
    vertices = vertices(1:3,:)';
    
    % Read faces
    fprintf('...Reading %d Faces\n',Nfaces);
    for f = 1:Nfaces,
        Nvert = fscanf(fid,'%d',1);
        faces(f,:) = fscanf(fid,'%d',Nvert);
    end
    % add 1 because GeomView vertices start at zero
    faces = faces + 1;
    
    fclose(fid);
    
    t = toc;
    fprintf('...done (%6.2f sec).\n',t);
    
end

return



% "dodec.off":
% 
% OFF
% 20 12 30
% 	1.214124 0.000000 1.589309
% 	0.375185 1.154701 1.589309
% 	-0.982247 0.713644 1.589309
% 	-0.982247 -0.713644 1.589309
% 	0.375185 -1.154701 1.589309
% 	1.964494 0.000000 0.375185
% 	0.607062 1.868345 0.375185
% 	-1.589309 1.154701 0.375185
% 	-1.589309 -1.154701 0.375185
% 	0.607062 -1.868345 0.375185
% 	1.589309 1.154701 -0.375185
% 	-0.607062 1.868345 -0.375185
% 	-1.964494 0.000000 -0.375185
% 	-0.607062 -1.868345 -0.375185
% 	1.589309 -1.154701 -0.375185
% 	0.982247 0.713644 -1.589309
% 	-0.375185 1.154701 -1.589309
% 	-1.214124 0.000000 -1.589309
% 	-0.375185 -1.154701 -1.589309
% 	0.982247 -0.713644 -1.589309
% 	5 0 1 2 3 4
% 	5 0 5 10 6 1
% 	5 1 6 11 7 2
% 	5 2 7 12 8 3
% 	5 3 8 13 9 4
% 	5 4 9 14 5 0
% 	5 15 10 5 14 19
% 	5 16 11 6 10 15
% 	5 17 12 7 11 16
% 	5 18 13 8 12 17
% 	5 19 14 9 13 18
% 	5 19 18 17 16 15
% 
% The "OFF" header tells us it's a polylist file. The second line in the
% file tells us that there are 20 vertices, 12 faces, and 30 edges. (The
% OOGL libraries presently don't use the edges value, so you can just use 0
% if you don't happen know the number of edges.) The next 20 lines give a
% list of vertices. The last 12 lines specify the faces: the first number
% is the number of vertices in that face. Since our polyhedron happens to
% be regular, all faces have the same number of vertices (in this case, 5).
% The rest of the numbers on the line are indices into the above list of
% vertices.

