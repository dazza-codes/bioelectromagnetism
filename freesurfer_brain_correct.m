function [FV_brain] = freesurfer_brain_correct(FV_brain,FV_pial)

% freesurfer_brain_correct - Correct intersections of brain/pial surfaces
%
% [FV_brain] = freesurfer_brain_correct(FV_brain,FV_pial)
%
% FV_brain.faces - Nx3
% FV_brain.vertices - Nx3
% FV_pial.faces - Nx3
% FV_pial.vertices - Nx3
%
% The FV_brain approximates an 'inner_skull', as it is the brain.tri
% file output by freesurfer during the skull strip procedure of the
% segmentation process.  This surface can have intersections
% with the pial surface, which is disturbing, so we should take
% some trouble to shift it outward along the surface normal to avoid
% intersections.  For practical purposes, we will have to do this
% regardless of the MRI intensity, so it is a detachment from the
% realistic geometry, unfortunatley.  We might assume that a 2 mm
% separation is a reasonable minimum.  The pial surface is assumed to
% be correct, so the 'brain' surface is moved outward along the surface
% normal, taking care not to introduce surface intersections in the brain
% surface itself (ie, constraining the inflation by neighbouring vertices).
% This is an iterative process for each vertex of the brain surface, each
% iteration moves a vertex and its neighbours by 0.5 mm.
%

% $Revision: 1.1 $ $Date: 2005/08/15 21:59:42 $

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

% History:  09/2005, Darren.Weber_at_radiology.ucsf.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $ $Date: 2005/08/15 21:59:42 $';
fprintf('FREESURFER_BRAIN_CORRECT [v %s]\n',ver(11:15));
tic

%fprintf('...sorry, under development!\n\n'); return

DIST = 2; % inner skull to pial surface separation is 2 mm

DIST_MOVE = 0.25; % move 0.25 mm

FV = FV_brain;

if isfield(FV,'edges'),
    if isempty(FV.edges),
        FV.edges = mesh_edges(FV);
    end
else
    FV.edges = mesh_edges(FV);
end


% get surface normals
fprintf('...estimating brain vertex surface normals\n');
hf = figure('Visible','off');
hp = patch('faces',FV.faces,'vertices',FV.vertices);
normals = get(hp,'VertexNormals');
close(hf); clear hf hp
%Convert to unit normals
[nrm,unit_normals] = colnorm(normals');
FV.normals = unit_normals';
clear nrm unit_normals

% % get surface normals
% fprintf('...calculating vertex surface normals: pial\n');
% hf = figure('Visible','off');
% hp = patch('faces',FV_pial.faces,'vertices',FV_pial.vertices);
% normals = get(hp,'VertexNormals');
% close(hf); clear hf hp
% %Convert to unit normals
% [nrm,normals] = colnorm(normals');
% FV_pial.normals = normals';
% clear nrm normals

indexSmallDist = ones(size(FV.vertices,1),1);

iteration = 0;
while indexSmallDist,

    iteration = iteration + 1;
    
    % find vertices in the reduced tesselation that match those of the dense
    % tesselation; this can take several hours to run with dsearchn!
    fprintf('...find pial vertices nearest to each brain vertex\n');
    indexPial4Brain = [];
    if exist('nearpoints','file'),
        [indexPial4Brain, bestDistSq] = nearpoints(FV.vertices',FV_pial.vertices');
        indexPial4Brain = indexPial4Brain';
        bestDist = sqrt(bestDistSq)';
    else
        [indexPial4Brain, bestDist] = dsearchn(FV_pial.vertices,FV.vertices);
        bestDist = sqrt( bestDist .* bestDist );
    end


    % find the brain vertices that are too close to the pial surface
    indexSmallDist = find(bestDist < DIST);
    
    % find all the immediate neighbours of these vertices
    move_indices = [];
    for i = 1:length(indexSmallDist),
        index = indexSmallDist(i);
        nearest = find(FV.edges(:,index));
        move_indices = [move_indices; index; nearest];
    end
    move_indices = unique(move_indices);

    % try to exclude the vertices of the cerebellum, 
    % by excluding vertex distances > 6 mm
    tmpDist = bestDist(move_indices);
    tmpIndex = find(tmpDist > (3 * DIST) );    
    move_indices = setdiff(move_indices, tmpIndex);
    clear tmpDist tmpIndex
    
    Nvert = length(move_indices);
    fprintf('...shifting %5d selected vertices outward %4.2f mm\n', Nvert, DIST_MOVE);

    % We have an irregular surface (ie, the inner skull surface) and it
    % may not have an origin at (0,0,0).  We have the surface vertex
    % normals, so now we move the vertex outward in the direction of
    % the vertex normal.  We simply add a multiple of the unit vertex
    % normal to the current vertex xyz location.
    
    move_vector = DIST_MOVE * FV.normals(move_indices,:);

    move_vertices = FV.vertices(move_indices,:);

    FV.vertices(move_indices,:) = move_vertices + move_vector;
    
    if iteration > 250,
        warning('cannot correct brain surface');
        return
    end
    
end

% apply one more move before smoothing
move_vector = DIST_MOVE * FV.normals;
move_vertices = FV.vertices;
FV.vertices = move_vertices + move_vector;

% smooth the surface, which can shrink it slightly
FV = mesh_vertex_smooth(FV);

FV_brain.vertices = FV.vertices;

t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code can be used to plot the vertices that are to be
% corrected
%
%     v = correct;
%     index = [scalpN oskullN];
%
%     patch('vertices',p.mesh.data.vertices{index(1)},'faces',p.mesh.data.faces{index(1)},...
%           'FaceColor',[.9 .9 .9],'Edgecolor','none','FaceAlpha',.6);
%     daspect([1 1 1]); axis tight; hold on
%
%     vert = p.mesh.data.vertices{index(1)};
%     x = vert(v,1);
%     y = vert(v,2);
%     z = vert(v,3);
%     plot3(x,y,z,'ro')
%
%     patch('vertices',p.mesh.data.vertices{index(2)},'faces',p.mesh.data.faces{index(2)},...
%           'FaceColor',[.0 .6 .0],'Edgecolor',[.2 .2 .2],'FaceAlpha',.8);
%
%     vert = p.mesh.data.vertices{index(2)};
%     x = vert(v,1);
%     y = vert(v,2);
%     z = vert(v,3);
%     plot3(x,y,z,'bo')
