function [FV] = mri_scalp_mesh(avw)

% MRI_SCALP_MESH: Find the scalp surface of an Analyze volume (avw)
%
% Useage: [FV] = mri_scalp_mesh(avw)
%
% FV is a struct with fields, vertices and faces.
%
%


FV = isosurface(avw.img,10);

return
