function reslice_axial(fname)
%
% function sag2ax(fname)
%
% Takes an analyze file acquired sagitally
% converted to SPM format using GE2SPM
% and rewrites it axially
%
% Souheil Inati
% Dartmouth College
% Jan. 2002
%

% Get the spm_vol structure
V = spm_vol(fname);

% Get the info in the header
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(fname);

% Make a new header with things swapped around
[dirname,imgname] = fileparts(V.fname);
ax_fname = fullfile(dirname,['ax_' imgname '.img']);
ax_DIM = [DIM(3) DIM(1) DIM(2)];
ax_VOX = [VOX(3) VOX(1) VOX(2)];
ax_SCALE = SCALE;
ax_TYPE = TYPE;
ax_OFFSET = OFFSET;
ax_DESCRIP = ['axialized ' DESCRIP];
ax_ORIGIN = [ORIGIN(3) DIM(1)-ORIGIN(1) DIM(2)-ORIGIN(2)];
spm_hwrite(ax_fname,ax_DIM,ax_VOX,ax_SCALE,ax_TYPE,ax_OFFSET,ax_ORIGIN,ax_DESCRIP);

V2 = spm_vol(ax_fname);

% Read in the data
[Y,XYZ] = spm_read_vols(V);

% Compute the orientation stuff
ROT = V.mat*inv(diag([VOX 1]));

% Reorient
data = permute(Y,[3 1 2]);
% Deal with GE image flip
data = flipdim(data,2);  % flip y
data = flipdim(data,3);  % flip z
% Deal with radiological flip
% data = flipdim(data,1);  % flip x  (make it neurological)

% Write out the reoriented data
spm_write_vol(V2,data);

return
