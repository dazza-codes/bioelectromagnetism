function sf_save_vrml(filename,points, edges, faces)
% SF_SAVE_VRML	Save surface as VRML file (for WebSpace)
% 		SF_SAVE_VRML(filename, points [, edges [, faces]]) 
%		saves the surface in the given file 
%		in as complete a form as possible

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES:
% - edges is a N x 2 matrix of point indices for endpoints
% - faces is a M x 3 matrix of point indices for corners

%%% THINGS TO DO:
% - check arguments for valid number and type
% - add VRML nodes for camera, lighting, material, text

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 2) help sf_save_vrml; return; end
if (nargin > 2)
    if (size(edges, 2) ~= 2)
	error(['Invalid number of columns in edge list (', ...
	       num2str(size(edges, 2)), ').']);
    end
end
if (nargin > 3)
    if (size(faces, 2) ~= 3)
	error(['Invalid number of columns in face list (', ...
	       num2str(size(faces,2)), ').']);
    end
end

%%% open file and write header
fp = fopen(filename,'w');
if (fp == -1)
    error(['Unable to open file ',filename,' for writing.']); 
end

fprintf(fp, '#VRML V1.0 ascii\n');
fprintf(fp, 'Separator {\n');

%%% coordinates of points (transformed array)
fprintf(fp, '  Coordinate3 { \n');
fprintf(fp, '    point [ \n');
fprintf(fp, '      %6.2f %6.2f %6.2f, \n', points');
fprintf(fp, '    ] \n');
fprintf(fp, '  } \n');

if (nargin == 2)
    %%% display points
    fprintf(fp, '  PointSet { \n');
    fprintf(fp, '  } \n');
    
elseif (nargin == 3)
    % convert to 0-based indexing
    edges = edges - 1;
    %%% display edges
    fprintf(fp, '  IndexedLineSet { \n');
    fprintf(fp, '    coordIndex [ \n');
    fprintf(fp, '      %d, %d, -1, \n', edges');
    fprintf(fp, '    ] \n');
    fprintf(fp, '  } \n');

elseif (nargin == 4)
    % convert to 0-based indexing
    faces = faces - 1;
    %%% display faces
    fprintf(fp, '  IndexedFaceSet { \n');
    fprintf(fp, '    coordIndex [ \n');
    fprintf(fp, '      %d, %d, %d, -1, \n', faces');
    fprintf(fp, '    ] \n');
    fprintf(fp, '  } \n');
    
else
    error(['Invalid number of input arguments (', num2str(nargin), ').']);
end

%%% write tail, close file, cleanup
fprintf(fp, '}\n');
fclose(fp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks
