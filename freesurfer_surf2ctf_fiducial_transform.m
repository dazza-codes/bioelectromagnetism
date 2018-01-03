
clear all
close all

if isunix,
  dataPath = '/data/freesurfer/subjects/ucsf_aw_fullcortex/';
else
  dataPath = 'e:\freesurfer\subjects\ucsf_aw_fullcortex\';
end

plot_fs_surf    = 0;
plot_fs_voxels  = 0;
plot_ctf_voxels = 0;

plot_fs_fid_ras = 0;
plot_ctf_fid    = 0;

plot_fs2ctf_fid  = 0;
plot_fs2ctf_surf = 0;


% the freesurfer MRI volume (cor-???) was converted to analyze
% format using the following command:
%  mri_convert orig/ -oid -1 0 0 -ojd 0 1 0 -okd 0 0 1 analyze/ucsf_aw_orig_axial_las.img

% This is what mri_convert calls the Analyze output matrix:
% avw_mat = [-1.000   0.000   0.000   129.000;
%             0.000   1.000   0.000  -129.000;
%             0.000   0.000   1.000  -129.000;
%             0.000   0.000   0.000     1.000;
%           ];

% Perhaps, by application of this matrix to the freesurfer surface,
% we obtain the freesurfer vertex locations in the space of the
% Analyze volume

%avwFile = [dataPath,'mri/analyze/ucsf_aw_orig_axial_las.hdr'];
%avw = avw_read(avwFile);


%-------------------------------------------------------
% Read CTF MRI volume

if isunix,
  ctfFile = [dataPath,'mri/analyze/ucsf_aw_orig_axial.mri'];
else
  ctfFile = [dataPath,'mri\analyze\ucsf_aw_orig_axial.mri']
end

% CTF MRI default volume index has an origin at the left, anterior, superior
% voxel, such that Sag slices increase from left to right, Cor
% slices increase from anterior to posterior, and Axi slices
% increase from superior to inferior

% CTF volume is 256^3 voxels, 1mm^3 each
% Sag increases from left to right (+X Right)
% Cor increases from anterior to posterior (+Y Posterior)
% Axi increases from superior to inferior (+Z Inferior)
% The origin is at the left, anterior, superior voxel
% CTF volume_index is: Sag, Cor, Axi (256^3, 1mm^3 voxels)

% The fiducial values can be read from the MRI volume, like this:
ctf.mri = ctf_read_mri(ctfFile);
headModel = ctf.mri.hdr.HeadModel_Info;

%-------------------------------------------------------
% Read CTF Head Shape

% if isunix,
%   ctfHeadFile = [dataPath,'mri/analyze/ucsf_aw_orig_axial.head.shape'];
% else
%   ctfHeadFile = [dataPath,'mri\analyze\ucsf_aw_orig_axial.head.shape']
% end
% HeadShape.vertices = ctf_read_headshape(ctfHeadFile);
% HeadShape.faces = convhulln(HeadShape.vertices);


%-------------------------------------------------------
% Read FreeSurfer Surface

if isunix,
  fsFile = [dataPath,'surf/rh.pial'];
else
  fsFile = [dataPath,'surf\rh.pial'];
end
[fsSurf.vertices,fsSurf.faces] = freesurfer_read_surf(fsFile);

% reduce the surface density for the sake of this testing
% fsSurf = reducepatch(fsSurf,10000);


%-------------------------------------------------------
% Convert FreeSurfer Surface into FreeSurfer voxel coordinates

fsVox.vertices = freesurfer_surf2voxels(fsSurf.vertices);

fsVox.faces = fsSurf.faces;

%-------------------------------------------------------
% Convert FreeSurfer voxels into CTF voxel coordinates

ctf.vertices.mri = freesurfer_surf2ctf(fsSurf.vertices);

ctf.faces.mri = fsSurf.faces;

%ctf_write_mrishape(ctf.vertices.mri,ctf.mri);

%-------------------------------------------------------
% Convert CTF voxels into CTF xyz coordinates

ctf.vertices.head = ctf_mri2head(ctf.vertices.mri,ctf.mri);

ctf.faces.head = fsSurf.faces;

%ctf_write_headshape(ctf.vertices.head,ctf.mri);


%-------------------------------------------------------
% Write out to BrainStorm

if isunix,
  dataPath = '/data/brainstorm_database/subjects/ucsf_aw/';
else
  dataPath = 'e:\brainstorm_database\subjects\ucsf_aw\';
end

subjecttessFile = [dataPath,'ucsf_aw_subjecttess'];


load( subjecttessFile )
plot_cortices = 0;
if plot_cortices,
    figure('name','current brainstorm')
    patch('vertices',Vertices{1}','faces',Faces{1},'facecolor',[.75 .7 .7],'edgecolor','none');
    view([0,0]); light
    figure('name','new brainstorm')
    patch('vertices',ctf.vertices.head / 100,'faces',ctf.faces.head,'facecolor',[.75 .7 .7],'edgecolor','none');
    view([0,0]); light
end

Comment = []; Comment{1} = 'cortex';
Faces = [];
Vertices = [];

Vertices{1} = ctf.vertices.head' / 100;  % convert cm to meter
Faces{1} = ctf.faces.head;

save(subjecttessFile, 'Comment', 'Faces', 'Vertices')

return











%-------------------------------------------------------
% Do some plots

if plot_fs_surf,
  
  % To confirm that the freesurfer surface is read correctly and the
  % fiducials are correct, we can plot them all
  figure('name','FreeSurfer Surface')
  hold on
  
  patch('faces',fsSurf.faces,'vertices',fsSurf.vertices, ...
    'edgecolor','none', ...
    'facecolor',[.8 .7 .7],'facealpha',0.5);
  
  scatter3(fs.volume_xyz.nas(1),fs.volume_xyz.nas(2),fs.volume_xyz.nas(3),40,'r','filled');
  scatter3(fs.volume_xyz.lpa(1),fs.volume_xyz.lpa(2),fs.volume_xyz.lpa(3),40,'g','filled');
  scatter3(fs.volume_xyz.rpa(1),fs.volume_xyz.rpa(2),fs.volume_xyz.rpa(3),40,'b','filled');
  text(fs.volume_xyz.nas(1),fs.volume_xyz.nas(2),fs.volume_xyz.nas(3),'. nasion');
  text(fs.volume_xyz.lpa(1),fs.volume_xyz.lpa(2),fs.volume_xyz.lpa(3),'. left');
  text(fs.volume_xyz.rpa(1),fs.volume_xyz.rpa(2),fs.volume_xyz.rpa(3),'. right');
  
  set(gca,'DataAspectRatio',[1 1 1]);
  set(gca,'Xlim',[-128 128],'Ylim',[-128 128],'Zlim',[-128 128])
  rotate3d
  view(2)
  light
end

if plot_fs_voxels,
  
  % To confirm that the freesurfer surface is read correctly and the
  % fiducials are correct, we can plot them all
  figure('name','FreeSurfer Voxels')
  hold on
  
  patch('faces',fsVox.faces,'vertices',fsVox.vertices, ...
    'edgecolor','none', ...
    'facecolor',[.8 .7 .7],'facealpha',0.5);
  
  scatter3(fs.volume_index.nas(1),fs.volume_index.nas(2),fs.volume_index.nas(3),40,'r','filled');
  scatter3(fs.volume_index.lpa(1),fs.volume_index.lpa(2),fs.volume_index.lpa(3),40,'g','filled');
  scatter3(fs.volume_index.rpa(1),fs.volume_index.rpa(2),fs.volume_index.rpa(3),40,'b','filled');
  text(fs.volume_index.nas(1),fs.volume_index.nas(2),fs.volume_index.nas(3),'. nasion');
  text(fs.volume_index.lpa(1),fs.volume_index.lpa(2),fs.volume_index.lpa(3),'. left');
  text(fs.volume_index.rpa(1),fs.volume_index.rpa(2),fs.volume_index.rpa(3),'. right');
  
  set(gca,'DataAspectRatio',[1 1 1]);
  set(gca,'Xlim',[0 256],'Ylim',[0 256],'Zlim',[0 256])
  rotate3d
  view([0,0])
  light
  
end

if plot_ctf_voxels,
  
  % Now we can plot the CTF fiducials, which should be about 90
  % degree rotations of the FreeSurfer fiducials
  
  figure('name','CTF Voxels')
  hold on
  
  patch('faces',ctf.shape.faces,'vertices',ctf.shape.vertices, ...
    'edgecolor','none', ...
    'facecolor',[.8 .7 .7],'facealpha',0.5);
  
  scatter3(ctf.volume_index.nas(1),ctf.volume_index.nas(2),ctf.volume_index.nas(3),40,'r','filled');
  scatter3(ctf.volume_index.lpa(1),ctf.volume_index.lpa(2),ctf.volume_index.lpa(3),40,'g','filled');
  scatter3(ctf.volume_index.rpa(1),ctf.volume_index.rpa(2),ctf.volume_index.rpa(3),40,'b','filled');
  text(ctf.volume_index.nas(1),ctf.volume_index.nas(2),ctf.volume_index.nas(3),'. nasion');
  text(ctf.volume_index.lpa(1),ctf.volume_index.lpa(2),ctf.volume_index.lpa(3),'. left');
  text(ctf.volume_index.rpa(1),ctf.volume_index.rpa(2),ctf.volume_index.rpa(3),'. right');
  
  set(gca,'DataAspectRatio',[1 1 1]);
  set(gca,'Xlim',[0 256],'Ylim',[0 256],'Zlim',[0 256])
  rotate3d
  view([0,-90])
  light
end
