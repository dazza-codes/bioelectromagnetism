function [outHeader] = GE_reorientHeader(inHeader, orient)
%
% outHeader = GE_reorientHeader(inHeader, orient)
%
% Fixes the analyze header file based on the acquistion
%
% Souheil Inati
% Dartmouth College
% May 2000
% souheil.inati@dartmouth.edu
%

% Initialize the output Header
outHeader = inHeader;

% Get the dimensions
nX = inHeader.dim(2);
nY = inHeader.dim(3);
nZ = inHeader.dim(4);
pX = inHeader.pixDim(2); % X Voxel Size
pY = inHeader.pixDim(3); % Y Voxel Size
pZ = inHeader.pixDim(4); % Z Voxel Size

% Reshape and flip around to match what SPM expects
% ie axial in radiological convention
% Orientation is 1=axial, 2=sagittal, 3=coronal
% with opposite sign if backwards slice order
switch orient
case  0 % Undefined Orientation
   fprintf('Orientation is undefined.  Check your images!');

case {1, -1} % Sagittal
   outHeader.dim(2) = nZ;
   outHeader.dim(3) = nX;
   outHeader.dim(4) = nY;
   outHeader.pixDim(2) = pZ;
   outHeader.pixDim(3) = pX;
   outHeader.pixDim(4) = pY;

case {2, -12}  % Axial
   outHeader.dim(2) = nX;
   outHeader.dim(3) = nY;
   outHeader.dim(4) = nZ;
   outHeader.pixDim(2) = pX;
   outHeader.pixDim(3) = pY;
   outHeader.pixDim(4) = pZ;

case {3, -3} % Coronal
   outHeader.dim(2) = nX;
   outHeader.dim(3) = nZ;
   outHeader.dim(4) = nY;
   outHeader.pixDim(2) = pX;
   outHeader.pixDim(3) = pZ;
   outHeader.pixDim(4) = pY;

end

% Set the origin to the center of the volume
outHeader.originator = [floor(outHeader.dim(2)/2) floor(outHeader.dim(3)/2) floor(outHeader.dim(4)/2) 0 0];

return
