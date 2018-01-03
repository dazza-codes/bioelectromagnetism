clear

dweber_paths

subject = fullfile(FS.subjectsDIR,'ucsf_cj','');

pial = freesurfer_read_subject(subject,{'pial','curv'});
pial.lh.pial.faces = pial.lh.pial.faces(:,[1 3 2]);
pial.rh.pial.faces = pial.rh.pial.faces(:,[1 3 2]);

surf = freesurfer_surf_combine(pial.lh.pial, pial.rh.pial);
surf.curv = [pial.lh.curv; pial.rh.curv];

reduceSurf = 0;

if reduceSurf,
  
  fprintf('...running reducepatch\n');
  surfReduced = reducepatch(surf, 80000);
  
  % find vertices in the reduced tesselation that match those of the dense
  % tesselation it's fast with nearpoints, but this can take several hours to
  % run with dsearchn!
  indexSparseInPial = [];
  if exist('nearpoints','file'),
    fprintf('...running nearpoints\n');
    indexSparseInPial = nearpoints(surfReduced.vertices',surf.vertices');
    indexSparseInPial = indexSparseInPial';
  else
    fprintf('...running dsearchn\n');
    indexSparseInPial = dsearchn(surf.vertices,surfReduced.vertices);
  end
  
  % assign the curvature from the dense tesselation
  % into the reduced tesselation
  surfReduced.curv = surf.curv(indexSparseInPial,:);

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

saveImages = 0;

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
  save_png('pialSurfL.png',Hf);
else
  pause(2)
end

view(viewR)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  save_png('pialSurfR.png',Hf);
else
  pause(2)
end

view(viewA)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  save_png('pialSurfA.png',Hf);
else
  pause(2)
end

view(viewP)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  save_png('pialSurfP.png',Hf);
else
  pause(2)
end

view(viewD)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  save_png('pialSurfD.png',Hf);
else
  pause(2)
end

view(viewV)
camlight(Hlight,'headlight');
drawnow
if saveImages,
  save_png('pialSurfV.png',Hf);
else
  pause(2)
end
