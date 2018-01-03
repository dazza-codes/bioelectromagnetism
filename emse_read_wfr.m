function [vertex,face,edge,meshtype,space] = emse_read_wfr(file,options)

% emse_read_wfr - read EMSE wireframe (.wfr)
%
% [vertex,face,edge,meshtype] = emse_read_wfr(file,[options])
%
% All values returned are in meters.  With the returned structures, 
% create a vertex & face matrix:
% 
%   vertex_matrix = [vertex.x; vertex.y; vertex.z]';
%   face_matrix = [face.vertex1; face.vertex2; face.vertex3]';
%
% These can be input to the patch command:
%
%    Hpatch = patch('Vertices',vertex_matrix,'Faces',face_matrix,...
%                   'EdgeColor',[.6 .6 .6],'FaceColor',[0.9 0.9 0.9]);
%
% See the patch and light commands to colour this object.
% 
% 'options' ... a cell array of strings.  By default it contains
% options = {'vertex','face','edge'}.  By default, this routine
% reads all available data from the emse file.  If 'options' is given
% and it doesn't contain one of these strings, that data will not be
% read or returned.
%
% meshtype is: 'unknown','scalp','outer skull','inner skull', or 'cortex'
%
% space - 'hspace' for head space (electrodes)
%         'vspace' for MRI volume space
%


% $Revision: 1.5 $ $Date: 2007/11/29 21:36:43 $

% Copyright (C) 2001  Darren L. Weber
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

% History:  10/98 Abbas Kouzani (egazk@flinders.edu.au)
%           09/01 Darren.Weber_at_radiology.ucsf.edu
%                 - created function, rather than script
%                 - added functionality to handle different
%                   minor revisions.
%           03/02 Darren.Weber_at_radiology.ucsf.edu
%                 - optimised fscanf processes for matlab, the
%                   function now operates in less than 1/2 time,
%                   especially for revision 3 data.
%           12/04 Darren.Weber_at_radiology.ucsf.edu
%                 - modified function name from mesh_emse2matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.5 $ $Date: 2007/11/29 21:36:43 $';
fprintf('\nEMSE_READ_WFR [v%s]\n',ver(11:15));

if ~exist('options','var'),
    options = {'vertex','face','edge'};
end

[path,name,ext] = fileparts(file);
file = fullfile(path,[name ext]);

[fid,msg] = fopen(file,'r');
if ~isempty(msg), error(msg); end

tic;

% Read prolog
version   =fscanf(fid,'%f',1);
file_type =fscanf(fid,'%f',1);
minor_rev =fscanf(fid,'%f',1);

fprintf('...WFR version = %d, minor_revision = %d, file-type = %d\n',...
    version,minor_rev,file_type);

if ~(file_type==8000 | file_type==4000)
    msg = sprintf('cannot read WFR file type: %d',file_type);
    error(msg);
end

% Read header (format depends on minor revision)
if(minor_rev==3)
  type      =fscanf(fid,'%f',1);
else
  radius    =fscanf(fid,'%f',1);
  vertex_no =fscanf(fid,'%d',1);
  face_no   =fscanf(fid,'%d',1);
  edge_no   =fscanf(fid,'%d',1);
  
  fprintf('...radius = %f meters\n',radius);
  fprintf('...reading %d vertices\n...%d faces\n...%d edges\n',...
          vertex_no,face_no,edge_no);

  if(minor_rev==1),
    type = 0;
  else,
    type = fscanf(fid,'%f',1);
  end
end

if(minor_rev==1)
  meshtype = 'unknown';
else
  if type >= 80000,
    space = 'vspace';                   % MRI "volume space"
    type = type - 80000;
  else
    space = 'hspace';                   % electrode "head space"
  end
  switch type
   case 0,
    meshtype = 'unknown';
   case { 64,  40},
    meshtype = 'scalp';
   case {128,  80},
    meshtype = 'outer skull';
   case {256, 100},
    meshtype = 'inner skull';
   case {512, 200},
    meshtype = 'cortex';
   otherwise,
    meshtype = 'unknown';
  end
end

fprintf('...mesh type: %s\n', meshtype);

% Read data (format depends on minor revision)
vertex = struct([]);
face = struct([]);
edge = struct([]);

if(minor_rev==3)
    
    % Read the whole file
    fprintf('...reading minor revision 3 data');
    Tmp = fscanf(fid,'%s%f%f%f',[4,inf]);
    fprintf('...done\n');
    
    % Vertices
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...creating vertex struct');
        % first get the numeric code for 'v'
        vcode = double('v');
        % now find all columns of Tmp with vcode in 1st row
        vindex = find(Tmp(1,:) == vcode);
        
        tmp = Tmp(2:4,vindex);
        tmp = num2cell(tmp);
        index = num2cell(1:size(tmp,2)); % start at 1 for Matlab index
        tmp = [index; tmp];
        vertex = cell2struct(tmp, {'index','x','y','z'});
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...skipping vertices\n');
    end
    
    % Faces
    if(max(strcmp('face',options)) > 0),
        fprintf('...creating face struct');
        % first get the numeric code for 't'
        tcode = double('t');
        % now find all columns of Tmp with tcode in 1st row
        tindex = find(Tmp(1,:) == tcode);
        
        % matlab vertex indices start at one,
        % not zero, so we add one to these emse values
        tmp = Tmp(2:4,tindex) + 1;
        tmp = num2cell(tmp);
        index = num2cell(1:size(tmp,2)); % start at 1 for Matlab index
        tmp = [index; tmp];
        face = cell2struct(tmp, {'index','vertex1','vertex2','vertex3'});
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...skipping faces\n');
    end
    
    % Edges
    fprintf('...no edges for minor revision 3\n');
    clear Tmp;
    
elseif(minor_rev~=4)
    
    % minor revision 1 & 2 format
    disp('...reading minor revision 1 or 2 data');
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...reading %d vertices',vertex_no);
        
        tmp = zeros(13,vertex_no);
        tmp = fscanf(fid,'%d%x%d%d%g%g%g%g%g%g%g%g%g',[13,vertex_no]);
        
        tmp(1,:) = tmp(1,:) + 1;        % add 1 for a Matlab index
        tmp      = num2cell(tmp);
        vertex   = struct(...
            'index',        tmp( 1,:),...
            'address',      tmp( 2,:),...
            'channel_index',tmp( 3,:),...
            'x',            tmp( 5,:),...
            'y',            tmp( 6,:),...
            'z',            tmp( 7,:),...
            'xnormal',      tmp( 9,:),...
            'ynormal',      tmp(10,:),...
            'znormal',      tmp(11,:),...
            'potential',    tmp(12,:),...
            'curvature',    tmp(13,:));
        clear tmp;
        fprintf('...done\n');
    else
        fprintf('...skipping vertices\n');
        tmp = zeros(13,vertex_no);
        tmp = fscanf(fid,'%d%x%d%d%g%g%g%g%g%g%g%g%g',[13,vertex_no]);
        clear tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('face',options)) > 0),
        fprintf('...reading %d faces',face_no);
        
        tmp = zeros(18,face_no);
        tmp = fscanf(fid,'%d%x%g%g%g%g%g%g%g%g%g%g%x%x%x%x%x%x',[18,face_no]);
        tmp(1,:) = tmp(1,:) + 1;        % start at 1 for Matlab index
        tmp      = num2cell(tmp);
        face     = struct(...
            'index',        tmp( 1,:),...
            'address',      tmp( 2,:),...
            'solid_angle',  tmp( 3,:),...
            'magnitude',    tmp( 4,:),...
            'potential',    tmp( 5,:),...
            'area',         tmp( 6,:),...
            'center_x',     tmp( 7,:),...
            'center_y',     tmp( 8,:),...
            'center_z',     tmp( 9,:),...
            'normal_x',     tmp(10,:),...
            'normal_y',     tmp(11,:),...
            'normal_z',     tmp(12,:),...
            'vertex1',      tmp(13,:),...
            'vertex2',      tmp(14,:),...
            'vertex3',      tmp(15,:),...
            'edge1',        tmp(16,:),...
            'edge2',        tmp(17,:),...
            'edge3',        tmp(18,:));
        clear tmp;
        fprintf('...done\n');
        
        % In minor rev4, the face vertex and edges
        % refer to the vertex.address field, so this
        % is corrected here.  Not sure how to avoid 'for'
        fprintf('...converting face vertices from address to index (this takes a while)');
        for i=1:face_no,
            face(i).vertex1 = find([vertex.address] == face(i).vertex1);
            face(i).vertex2 = find([vertex.address] == face(i).vertex2);
            face(i).vertex3 = find([vertex.address] == face(i).vertex3);
        end
        fprintf('...done\n');
    else
        disp('...skipping faces');
        tmp = zeros(18,face_no);
        tmp = fscanf(fid,'%d%x%g%g%g%g%g%g%g%g%g%g%x%x%x%x%x%x',[18,face_no]);
        clear tmp;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('edge',options)) > 0),
        fprintf('...reading %d edges',edge_no);
        
        tmp = zeros(4,edge_no);
        tmp = fscanf(fid,'%d%x%x%x',[4,edge_no]);
        
        tmp(1,:) = tmp(1,:) + 1;        % start at 1 for Matlab index
        tmp      = num2cell(tmp);
        edge     = struct(...
            'index',        tmp(1,:),...
            'address',      tmp(2,:),...
            'vertex1',      tmp(3,:),...
            'vertex2',      tmp(4,:));
        clear tmp;
        fprintf('...done\n');
        fprintf('...converting edge vertices from address to index (this takes a while)');
        for i=1:edge_no,
            edge(i).vertex1 = find([vertex.address] == edge(i).vertex1);
            edge(i).vertex2 = find([vertex.address] == edge(i).vertex2);
        end
        fprintf('...done\n');
        
        if(max(strcmp('edge',options)) > 0),
            fprintf('...converting face edges from address to index (this takes a while)');
            for i=1:face_no,
                face(i).edge1 = find([edge.address] == face(i).edge1);
                face(i).edge2 = find([edge.address] == face(i).edge2);
                face(i).edge3 = find([edge.address] == face(i).edge3);
            end
            fprintf('...done\n');
        end
    else
        disp('...skipping edges');
        %for i=1:edge_no,
        %    tmp=fscanf(fid,'%*d',1);
        %    tmp=fscanf(fid,'%*x',3);
        %end
    end
else
    % minor revision 4 format
    disp('...Reading Minor Revision 4 Data');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('vertex',options)) > 0),
        fprintf('...reading %d vertices',vertex_no);
        
        % start at 1 for Matlab index
        index = meshgrid(1:1:vertex_no,1);
        index = num2cell(index);
        
        tmp = zeros(11,vertex_no);
        tmp = fscanf(fid,'%d%d%g%g%g%d%g%g%g%g%g',[11,vertex_no]);
        
        tmp      = num2cell(tmp);
        vertex   = struct(...
            'index',        index,...
            'channel_index',tmp( 1,:),...
            'x',            tmp( 3,:),...
            'y',            tmp( 4,:),...
            'z',            tmp( 5,:),...
            'xnormal',      tmp( 7,:),...
            'ynormal',      tmp( 8,:),...
            'znormal',      tmp( 9,:),...
            'potential',    tmp(10,:),...
            'curvature',    tmp(11,:));
        clear tmp index;
        fprintf('...done\n');
    else
        disp('...skipping vertices');
        tmp = zeros(11,vertex_no);
        tmp = fscanf(fid,'%d%d%g%g%g%d%g%g%g%g%g',[11,vertex_no]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('face',options)) > 0),
        fprintf('...reading %d faces',face_no);
        
        % start at 1 for Matlab index
        index = meshgrid(1:1:face_no,1);
        index = num2cell(index);
        
        tmp = zeros(16,face_no);
        tmp = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%d%d%d%d%d%d',[16,face_no]);
        
        tmp(11:16,:) = tmp(11:16,:) + 1;
        
        tmp      = num2cell(tmp);
        face     = struct(...
            'index',        index,...
            'solid_angle',  tmp( 1,:),...
            'magnitude',    tmp( 2,:),...
            'potential',    tmp( 3,:),...
            'area',         tmp( 4,:),...
            'center_x',     tmp( 5,:),...
            'center_y',     tmp( 6,:),...
            'center_z',     tmp( 7,:),...
            'normal_x',     tmp( 8,:),...
            'normal_y',     tmp( 9,:),...
            'normal_z',     tmp(10,:),...
            'vertex1',      tmp(11,:),...
            'vertex2',      tmp(12,:),...
            'vertex3',      tmp(13,:),...
            'edge1',        tmp(14,:),...
            'edge2',        tmp(15,:),...
            'edge3',        tmp(16,:));
        clear tmp index;
        fprintf('...done\n');
    else
        disp('...skipping faces');
        tmp = zeros(16,face_no);
        tmp = fscanf(fid,'%g%g%g%g%g%g%g%g%g%g%d%d%d%d%d%d',[16,face_no]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(max(strcmp('edge',options)) > 0),
        fprintf('...reading %d edges',edge_no);
        
        % start at 1 for Matlab index
        index = meshgrid(1:1:edge_no,1);
        address = num2cell(index * 0);
        index = num2cell(index);
        
        tmp = zeros(2,edge_no);
        tmp = fscanf(fid,'%d%d',[2,edge_no]);
        tmp = tmp + 1;
        tmp = num2cell(tmp);
        
        edge = struct(...
            'index',   index,...
            'address', address,...
            'vertex1', tmp(1,:),...
            'vertex2', tmp(2,:));
        clear tmp index address;
        fprintf('...done\n');
    else
        fprintf('...skipping edges\n');
        %tmp=fscanf(fid,'%*d',2*edge_no);
    end
end
fclose(fid);

t=toc; fprintf('...done (%5.2f sec)\n\n',t);

return
