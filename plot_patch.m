function plot_patch(FV,varargin);
%PLOT_PATCH - plots a vertex or face and its neighbors
% function plot_patch(FV,varargin);
% The arguments should be entered in pairs of:
%   'faces',indices,
%   'vertices',indices,
%   'figure', handle,
%   'all', ignore
%
% If 'faces', then each face in the index is painted green with its immediate
% face neighbors.
% If 'vertices', then each vertex in the index is marked green with its
%   immediate face neighbors.
% if 'figure', then the figure handle is used, otherwise a new one is made.
% If 'all', then the entire FV is displayed as a semitransparent gray
%   background, and indices are ignored
% 
% See also TESSELLATION_STATS

%<autobegin> ---------------------- 21-May-2004 15:14:19 -----------------------
% --------- Automatically Generated Comments Block Using AUTO_COMMENTS ---------
%
% CATEGORY: Visualization
%
% At Check-in: $Author: psdlw $  $Revision: 1.1 $  $Date: 2004/11/12 01:32:35 $
%
% This software is part of BrainStorm Toolbox Version 2.0 (Alpha) 19-May-2004
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2004 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 21-May-2004 15:14:19 -----------------------

% ----------------------------- Script History ---------------------------------
% JCM 15-May-2004 Creation
% ----------------------------- Script History ---------------------------------



METHODS = varargin(1:2:end);
Indices = varargin(2:2:end);


% first calculate what faces go with which vertex

numTri = size(FV.faces,1); % number of faces
numVert = size(FV.vertices,1); % number of vertices

[VertexNumbering,I] = sort(FV.faces(:)); % sorted Vertex numbers

FacesNumbering = rem(I-1,numTri)+1; % triangle number for each Vertex

% For the ith vertex, then FacesNumbering(VertexNumbering == i) returns the indices of the
%  polygons attached to the ith vertex.
%
% For the set of vertices in the row vector sv (e.g sv = [3 5 115 121]), then use
%  [i,ignore] = find(VertexNumbering(:,ones(1,length(sv))) == (ones(size(VertexNumbering,1),1)*sv));
%  (compares the Vertex numbers to the indices, extracts the row indices into i)
%  then FacesNumbering(i) returns the indices. Apply unique to clean up.

% So now we know what faces are connected to each vertex


% --------------------------- Triangle Statistics ------------------------------
% calculate the area and normals for each triangle

Vertices = FV.vertices(FV.faces',:);
% each set of three rows of Vertices is one triangle

% now difference them to get the vectors on two sides
dVertices = diff(Vertices);
dVertices(3:3:end,:) = []; % remove the transition between triangles

% now each pair of rows in dVertices represents each triangle
% row 1 is vector from 1 to 2
% row 2 is vector from 2 to 3

% v1 = dVertices(1:2:end,:)'; % each column is side one of a triangle
% v2 = dVertices(2:2:end,:)'; % side two

% right-hand rule, counter-clockwise ordering of the triangle yields a positive
% upward area and normal.
% Call a fast subfunction of this function:
WeightedNormals = cross(dVertices(1:2:end,:)',dVertices(2:2:end,:)')/2; 
% each column is the normal for each triangle
% the length the vector gives the area

FaceArea = sqrt(sum(WeightedNormals .* WeightedNormals)); % the area
FaceNormal = WeightedNormals ./ (FaceArea([1 1 1],:)); % normalize them

% now calculate the centers of each triangle
FaceCenter = cumsum(Vertices);
FaceCenter = FaceCenter(3:3:end,:); % every third summation for every triangle
FaceCenter = diff([0 0 0;FaceCenter])'/3; % average of each summation
% each column is the mean of the vertices of the triangles
% so now we know the center, area, and the normal vector of each triangle


ih = find(strcmp('figure',varargin));
if ~isempty(ih),
   hf = figure(varargin{ih+1});
else
   hf = figure; % open new figure
end

hold on

for imethod = 1:length(METHODS),
   
   switch lower(METHODS{imethod})
      case 'faces'
         
         iFace = FV.faces([Indices{imethod}],:);
         iFace = iFace(:)'; % ensure row vector
         
         % get the faces attached to these vertices
         % For the set of vertices in the row vector sv (e.g sv = [3 5 115 121]), then
         % use
         %  [i,ignore] = find(VertexNumbering(:,ones(1,length(sv))) == (ones(size(VertexNumbering,1),1)*sv));
         %  (compares the Vertex numbers to the indices, extracts the row indices into i)
         %  then FacesNumbering(i) returns the indices.
         
         [fndx, ignore] = find(VertexNumbering(:,ones(1,length(iFace))) == (ones(size(VertexNumbering,1),1)*iFace));
         fndx = FacesNumbering(fndx); % the faces attached to these vertices
         fndx = unique(fndx); % ensure unique
         nifndx = fndx(:)';
         nifndx(intersect(nifndx,iFace)) = []; % remove the ith face
         
         % not the ith face
         h = patch('vertices',FV.vertices,'faces',FV.faces(nifndx,:),'facecolor','r','edgecolor','k');
         % the ith face
         hi = patch('vertices',FV.vertices,'faces',FV.faces([Indices{imethod}],:),'facecolor','g');
         plot3(FaceCenter(1,fndx),FaceCenter(2,fndx),FaceCenter(3,fndx),'*')
         quiver3(FaceCenter(1,fndx),FaceCenter(2,fndx),FaceCenter(3,fndx),...
            FaceNormal(1,fndx),FaceNormal(2,fndx),FaceNormal(3,fndx),0.25);
         set(h,'FaceAlpha',.8)
         
      case 'vertices'
         % plot each vertex with a green point and the immediate faces around it
         
         % get the vertices
         iVert = [Indices{imethod}];
         [fndx, ignore] = find(VertexNumbering(:,ones(1,length(iVert))) == (ones(size(VertexNumbering,1),1)*iVert));
         fndx = FacesNumbering(fndx); % the faces attached to these vertices
         fndx = unique(fndx); % ensure unique
         
         h = patch('vertices',FV.vertices,'faces',FV.faces(fndx,:),'facecolor','r','edgecolor','k');
         hold on
         plot3(FaceCenter(1,fndx),FaceCenter(2,fndx),FaceCenter(3,fndx),'*');
         plot3(FV.vertices(iVert,1),FV.vertices(iVert,2),FV.vertices(iVert,3),'g+');
         ma = mean(FaceArea(fndx)); % mean area
         quiver3(FaceCenter(1,fndx),FaceCenter(2,fndx),FaceCenter(3,fndx),...
            FaceNormal(1,fndx),FaceNormal(2,fndx),FaceNormal(3,fndx),.25);
         set(h,'FaceAlpha',.8)
         
      case 'all'
         h = patch(FV,...
            'facecolor', [.9 .9 .9],'edgecolor',[.8 .8 .8]);
         set(h,'facealpha',.5,'edgealpha',.5);
      otherwise
         % unknown method
   end
end

axis equal
axis vis3d

hold off
cameratoolbar('Show'); % activate the camera toolbar
ret = cameratoolbar; % for strange reason, this activates the default orbit mode.
drawnow
