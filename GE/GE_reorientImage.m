function [outVol] = GE_reorientImage(inVol, orient)
%
%GE_reorientImage
%
% function [outVol] = GE_reorientImage(inVol, orient)
% reorients the inVol 3D volume to be SPM compatible
% based on the orient flag
%
% Souheil J. Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

[nX nY nZ] = size(inVol);

% Reshape and flip around to match what SPM expects
% ie axial in radiological convention
% Orientation is 1=axial, 2=sagittal, 3=coronal
% with opposite sign if backwards slice order
%
% This has been tested with Ax,Sag,Cor with slices going
% both ways.  Also for Oblique axial.
% Don't count on double obliques or anything really fancy
%


switch orient
case  0 % Undefined Orientation
   fprintf('Orientation is undefined.  Check your images!');
   outArray = inArray;

case  1 % Sagittal (L to R)
   outVol = permute(inVol,[3 1 2]);
   outVol(1:nZ,:,:) = outVol(nZ:-1:1,:,:);
   outVol(:,1:nX,:) = outVol(:,nX:-1:1,:);
   outVol(:,:,1:nY) = outVol(:,:,nY:-1:1);

case -1 % Sagittal (R to L)
   outVol = permute(inVol,[3 1 2]);
   outVol(:,1:nX,:) = outVol(:,nX:-1:1,:);
   outVol(:,:,1:nY) = outVol(:,:,nY:-1:1);

case  2  % Axial (I to S)
   outVol = inVol;
   outVol(:,1:nY,:) = outVol(:,nY:-1:1,:);

case -2  % Axial (S to I) 
   outVol = inVol;
   outVol(:,1:nY,:) = outVol(:,nY:-1:1,:);
   outVol(:,:,1:nZ) = outVol(:,:,nZ:-1:1);

case  3 % Coronal (A to P)
   outVol = permute(inVol,[1 3 2]);
   outVol(:,:,1:nY) = outVol(:,:,nY:-1:1);
   outVol(:,1:nZ,:) = outVol(:,nZ:-1:1,:);

case -3 % Coronal (P to A)
   outVol = permute(inVol,[1 3 2]);
   outVol(:,:,1:nY) = outVol(:,:,nY:-1:1);

end

return
