function freesurfer_surf_view(subjectPath,subjectID,surfName,saveImages)

% freesurfer_surf_view - view a freesurfer surface
%
% freesurfer_surf_view(subjectPath,subjectID,surfName,saveImages)
%
% subjectPath - full path to SUBJECTS_DIR
% subjectID - subject to view (e.g., 'average7')
% surfName - surface name to view (e.g., 'pial')
% saveImages - save .png snapshots of each view (default = 0)
%              images are saved into [surfName,'_view?.png']
%              within subjectPath/subjectID/tiff/
%

% $Revision: 1.1 $ $Date: 2007/03/09 22:24:25 $

% Copyright (C) 2007  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation; either version 2 of the License, or (at your option)
% any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% for more details.
% 
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc., 
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

% History:  03/2007, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ver = '$Revision: 1.1 $ $Date: 2007/03/09 22:24:25 $';
fprintf('FREESURFER_SURF_VIEW [v %s]\n',ver(11:15));

if(nargin < 3)
  help freesurfer_surf_view;
  error('...too few arguments')
end

if ~exist('saveImages)','var'),
  saveImages = 0;
end
if isempty(saveImages),
  saveImages = 0;
end

subjectPath = fullfile(subjectPath,subjectID,'');
imagePath = fullfile(subjectPath,'tiff','');

fs = freesurfer_read_subject(subject,{surfName,'curv','thickness'});
surf = freesurfer_surf_combine(fs.lh.(surfName), fs.rh.(surfName));
surf.faces = surf.faces(:,[1 3 2]);
surf.curv = [fs.lh.curv; fs.rh.curv];
clear fs

reduceSurf = 0;
if reduceSurf,
  
  fprintf('...running reducepatch\n');
  surfReduced = reducepatch(surf, 80000);
  
  % find vertices in the reduced tesselation that match those of the dense
  % tesselation it's fast with nearpoints, but this can take several hours to
  % run with dsearchn!
  indexSparseInFull = [];
  if exist('nearpoints','file'),
    fprintf('...running nearpoints\n');
    indexSparseInFull = nearpoints(surfReduced.vertices',surf.vertices');
    indexSparseInFull = indexSparseInFull';
  else
    fprintf('...running dsearchn\n');
    indexSparseInFull = dsearchn(surf.vertices,surfReduced.vertices);
  end
  
  % assign the curvature from the dense tesselation
  % into the reduced tesselation
  surfReduced.curv = surf.curv(indexSparseInFull,:);

  [Hf,Hp] = freesurfer_plot_curv(surfReduced, surfReduced.curv)
  
else
  
  [Hf,Hp] = freesurfer_plot_curv(surf, surf.curv)
  
end

%surf = mesh_smooth_vertex(surf);
%% How do we recalculate the surface curvature?
%[Hf,Hp] = freesurfer_plot_surf([],surf)

colormap(flipud(gray(200)))


colorbar off


% plot view config:
daspect([1,1,1]);
camproj perspective 
camva(7)
camtarget([0,0,0])

Hlight = camlight('headlight');

%lighting phong
%set(gcf,'Renderer','zbuffer')
%lighting gouraud
%set(gcf,'Renderer','OpenGL')

drawnow


viewL = [-90,  0];
viewR = [ 90,  0];
viewA = [180,  0];
viewP = [  0,  0];
viewD = [  0, 90];
viewV = [  0,-90];

view(viewL)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewL.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end

view(viewR)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewR.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end

view(viewA)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewA.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end

view(viewP)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewP.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end

view(viewD)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewD.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end

view(viewV)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  imageFile = fullfile(imagePath,[surf,'_viewV.png']);
  save_png(imageFile,Hf);
else
  pause(2)
end
