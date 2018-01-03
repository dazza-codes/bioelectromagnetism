function [pvc]=PVC_eqW(b,R,PSFx,PSFy,PSFz)

% Partial Volume Correction 
% 
% Usage : [pvc]=PVC(b,R,PSFx,PSFy,PSFz)
%
% Inputs
% b    = Measured PET Data (3d matrix) b(1:xdim,1:ydim,1:zdim) (b must be in counts/voxel)
% R    = ROI Definitions (4d binary matrix) R(1:xdim,1:ydim,1:zdim,1:no_regions)
% PSFx = Matrix Point Spread Function in x Direction. The ith row is the kernel for the ith voxel.
% PSFy = Matrix Point Spread Function in y Direction. The ith row is the kernel for the ith voxel.
% PSFz = Matrix Point Spread Function in z Direction. The ith row is the kernel for the ith voxel.
%
% Outputs
% These are contained in the structure pvc
% pvc.xm is the mean parameter estimate for the true activity concentrations in the ROIs
% pvc.xv is contains the associated variance estimates  
% pvc.ROIsize is size of each region.
%
% Tensor Implementation of the weighted least squares Partial Volume Correction Problem.
% PSF values are for reconstructed voxel size of b. Each number represents the recovery 
% coefficient in the next voxel for a point source in position 1 of the vector. 
%
% This method applies equal weighting to each voxel
% 
% (c) Roger Gunn & John Aston, 26-02-2000


pvc.xm=[];pvc.xv=[];
if nargin~=5 disp('Usage : [xm,xv]=PVC(b,R,PSFx,PSFy,PSFz)');return;end
if length(size(R))~=4;disp('Error: ROI matrix should be 4D');return;end
if length(size(b))~=3;disp('Error: b matrix should be 3D');return;end
if size(R,4)<2;disp('Error: There must be at least 1 ROI');return;end
if sum(size(R)==[size(b) 0])~=3;disp('Error: ROI and Image volume have different dimensions');return;end
if sum(size(PSFx)~=size(b,1));disp('Error: Check PSFx');return;end
if sum(size(PSFy)~=size(b,2));disp('Error: Check PSFy');return;end
if sum(size(PSFz)~=size(b,3));disp('Error: Check PSFz');return;end
if (size(b,1)==1|size(b,2)==1|size(b,3)==1);disp('Try a 3D Volume, It''s much nicer !');return;end

[xdim ydim zdim] 	= size(b);
no_regions 			= size(R,4);
PSFz 					= PSFz';


% Blur ROI's
for r=1:no_regions
   fprintf('.');
	for x=1:size(b,3);R(:,:,x,r)=PSFx(1:xdim,1:xdim)*squeeze(R(:,:,x,r));end
	for y=1:size(b,1);R(y,:,:,r)=PSFy(1:ydim,1:ydim)*squeeze(R(y,:,:,r));end
   for z=1:size(b,2);R(:,z,:,r)=squeeze(R(:,z,:,r))*PSFz(1:zdim,1:zdim);end
end

% xm = inv(R'P'PR)R'Pb
PR 		= reshape(R,xdim*ydim*zdim,no_regions);
pvc.xm	= inv(PR'*PR)*PR'*reshape(b,xdim*ydim*zdim,1);
pvc.xv	= inv(PR'*PR)*sum((reshape(b,xdim*ydim*zdim,1)-PR*pvc.xm).^2)/(xdim*ydim*zdim-no_regions);

pvc.xm=pvc.xm';
return;
