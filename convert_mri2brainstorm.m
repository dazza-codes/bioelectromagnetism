
% Script to convert freesurfer analyze files into brainstorm format
% The Freesurfer analyze files were created using 
% mri_convert -oid 1 0 0 -ojd 0 1 0 -okd 0 0 1 orig subject_orig_axial_ras.img
% to create an axial volume in with neurological orientation,
% ie, +x is Right, +y is Anterior, +z is Superior (RAS).

clear all

coregister = 0;
elecplot = 0;


% Fiducial points in RAS volume, obtained using avw_view (in meters)
% 
%                  Nasion                 Right                Left
%
mriFID.sub{ 1} = 'c01';
mriFID.xyz{ 1} = [  0.005  0.097 -0.022;  0.071  0.023 -0.052; -0.072  0.016 -0.050 ];
mriFID.sub{ 2} = 'c02';
mriFID.xyz{ 2} = [  0.007  0.088 -0.007;  0.081  0.015 -0.046; -0.070  0.012 -0.051 ];
mriFID.sub{ 3} = 'c03';
mriFID.xyz{ 3} = [ -0.002  0.098 -0.008;  0.055  0.040 -0.038; -0.095  0.030 -0.046 ];
mriFID.sub{ 4} = 'c04';
mriFID.xyz{ 4} = [  0.000  0.092  0.000;  0.070  0.024 -0.036; -0.076  0.031 -0.036 ];
mriFID.sub{ 5} = 'c05';
mriFID.xyz{ 5} = [  0.010  0.094 -0.010;  0.075  0.024 -0.048; -0.072  0.024 -0.052 ];
mriFID.sub{ 6} = 'c06';
mriFID.xyz{ 6} = [ -0.002  0.088  0.000;  0.077  0.008 -0.048; -0.075 -0.008 -0.048 ];
mriFID.sub{ 7} = 'c07';
mriFID.xyz{ 7} = [ -0.001  0.097 -0.018;  0.078  0.008 -0.055; -0.076  0.008 -0.055 ];
mriFID.sub{ 8} = 'c08';
mriFID.xyz{ 8} = [  0.000  0.096 -0.019;  0.087  0.005 -0.039; -0.065  0.004 -0.038 ];
mriFID.sub{ 9} = 'c09';
mriFID.xyz{ 9} = [  0.010  0.086 -0.020;  0.100 -0.020 -0.080; -0.060  0.020 -0.058 ];
mriFID.sub{10} = 'c10';
mriFID.xyz{10} = [  0.000  0.091 -0.018;  0.076  0.011 -0.052; -0.068  0.006 -0.052 ];

mriFID.sub{11} = 'p02';
mriFID.xyz{11} = [  0.002  0.091  0.016;  0.080  0.025 -0.036; -0.080  0.015 -0.043 ];
mriFID.sub{12} = 'p04';
mriFID.xyz{12} = [ -0.001  0.086 -0.009;  0.086  0.002 -0.062; -0.070  0.003 -0.065 ];
mriFID.sub{13} = 'p05';
mriFID.xyz{13} = [ -0.002  0.094 -0.001;  0.068  0.025 -0.040; -0.071  0.005 -0.042 ];
mriFID.sub{14} = 'p06';
mriFID.xyz{14} = [ -0.002  0.084  0.000;  0.082  0.013 -0.050; -0.066  0.013 -0.052 ];
mriFID.sub{15} = 'p07';
mriFID.xyz{15} = [  0.001  0.092  0.015;  0.080  0.003 -0.033; -0.070  0.004 -0.032 ];
mriFID.sub{16} = 'p08';
mriFID.xyz{16} = [ -0.003  0.095 -0.002;  0.070  0.018 -0.040; -0.074  0.022 -0.035 ];
mriFID.sub{17} = 'p09';
mriFID.xyz{17} = [ -0.002  0.100  0.004;  0.100  0.002 -0.028; -0.050  0.000 -0.036 ];





data = 'd:\matlab\brainstorm_v1\subjects\';

cd(data);

% Load data
for s = {'c01','c02','c03','c04','c05','c06','c07','c08','c09','c10',...
               'p02',      'p04','p05','p06','p07','p08','p09'},
    
       sub = sprintf('%s',char(s));
       cd(sub)
       
       p = eeg_toolbox_defaults;
       
       % Load the MRI volume (256^3, 1mm^3)
       p.mri.path = sprintf('d:\\freesurfer\\subjects\\ptsdpet-%s\\mri\\analyze\\',char(sub));
       p.mri.file = sprintf('%s_orig_axial_ras.img',char(sub));
       p.mri.plot = 0;
       p = mri_open(p);
       
       % -- Create the Patient Coordinate System (PCS) struct
       
       % Get surface fiducials from mriFID struct above
       Nfid = strmatch(sub,mriFID.sub);
       % 3x4, nasion, right, left, origin fiducial points in rows
       % Multiply the FID points by 1000 to get mm, rather than meters
       % Also add the origin (128)
       p.mriFID = [mriFID.xyz{Nfid} .* 1000 + 128; 128 128 128]';
       
       PCS.R = eye(3);     % [3x3 double] rotations
       PCS.t = zeros(3,1); % [3x1 double] translations
       PCS.Comment = 'NEUROMAG';
       PCS.PCSFiducial(:,1)  = p.mriFID(:,1); % NAS
       PCS.PCSFiducial(:,2)  = p.mriFID(:,3); % LPA
       PCS.PCSFiducial(:,3)  = p.mriFID(:,2); % RPA
       PCS.PCSFiducial(:,4)  = p.mriFID(:,4); % Origin
       PCS.CubeFiducial      = PCS.PCSFiducial;
       PCS.FiducialName = {'NAS'  'LPA'  'RPA'  'Origin'};
       
       % Load the scalp mesh
       p.mesh.path = sprintf('d:\\data_source\\%s\\meshes\\',char(sub));
       p.mesh.file = sprintf('%s_scalp.wfr',char(sub));
       p.mesh.type = 'emse';
       p = mesh_open(p);
       
       Scalp = p.mesh.data.vertices{1}';
       
       % The electrodes have already been coregistered and the transformed
       % coordinates saved, otherwise this code might do the job
       
       if coregister,
           
           % Get surface fiducials from mriFID struct above
           Nfid = strmatch(sub,mriFID.sub);
           % 3x3, nasion, right, left fiducial points in rows
           p.mriFID = mriFID.xyz{Nfid};
           
           % Load the electrode data
           p.elec.path = sprintf('d:\\data_source\\%s\\meshes\\',char(sub));
           p.elec.file = sprintf('%s_124fit.txt',char(sub));
           p.elec.plot = 0;
           p = elec_open(p);
           
           % Get electrode fiducials
           Efid = [p.elec.data.nasion; p.elec.data.rpa; p.elec.data.lpa];
           
           % Calculate coregistration transform
           T = elec_coregister(Efid,p.mriFID);
           % Create the PCS struct
           PCS.R = T([1:3;1:3]);
           PCS.t = T(4,1:3);
       end
       
       if elecplot,
           patch('vertices',p.mesh.data.vertices{1},'faces',p.mesh.data.faces{1},...
               'FaceColor',[1 0 0],'Edgecolor','none','FaceAlpha',.6);
           lighting phong, material dull, camlight headlight, hold on
           % Load the electrode data
           p.elec.path = sprintf('d:\\matlab\\brainstorm_v1\\studies\\%s\\',char(sub));
           p.elec.file = sprintf('%s_channel.mat',char(sub));
           p.elec.type = 'brainstorm';
           p.elec.plot = 0;
           p = elec_open(p);
           plot3(p.elec.data.x,...
                 p.elec.data.y,...
                 p.elec.data.z,'bo')
           daspect([1 1 1]); axis tight
           mouse_rotate
           return
       end
       
       
       % Maybe use FSL FAST result here?
       Segment = [];
       
       p.mri.data.fileprefix = sprintf('%s',char(sub));
       avw2brainstorm(p.mri.data,Segment,Scalp,PCS);
       
       cd ..
       
end
