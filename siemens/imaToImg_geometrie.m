function err = imaToImg_geometrie( imgFilename, path, params);
%
% Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
% calculate the rotations / translations and fix them in a .mat file
%
% intro / remember:
% coord in scanner: x == horizontal, to left is +
%                   y == vertical, up is +, therefor: z == to the feets is +
%
% coord in spm: rot about 180° versus y.
%
dir = pwd; cd( path); err=0;

M           = eye( 4);

Tspm        = eye( 3);  % transformation from MRT Koord to SPM world. 
Tspm(1,1)   = -1;   
Tspm(3,3)   = -1;   

Zoom  = zeros( 3);      % scaling in row, col and norm direction (Units == [mm]).
Zoom( 1, 1) = params.file.FoV(1) / params.file.matrix(1);   
Zoom( 2, 2) = params.file.FoV(2) / params.file.matrix(2);  
Zoom( 3, 3) = params.file.sliceThickness*(1+params.file.distFactor);

Tmrt     = zeros( 3);   % tranformation from native to scanned plot-orientation in MRT-KoordSys
Tmrt(1:3,1) = params.file.rowVect;  
Tmrt(1:3,2) = params.file.colVect;  
Tmrt(1:3,3) = params.file.normVect; 

M(1:3,1:3) = Tspm * Tmrt * Zoom; % total transformation 

% vector in MRT coord-sys to the voxel 1 1 1.
MRT_BBorigin = params.file.centerPoint - Zoom( 3, 3)/2*params.file.normVect - ...
   (params.file.matrix(1)+1)/2*Zoom( 1, 1) * params.file.rowVect - ...
   (params.file.matrix(2)+1)/2*Zoom( 2, 2) * params.file.colVect;

MRTcenter = params.file.centerPoint + Zoom( 3, 3) * (params.file.nSlices-1)/2 * params.file.normVect;

% this vector tells SPM where to locate the Voxel 1 1 1 (Units == [mm]).
M(1:3,4)    = Tspm * MRT_BBorigin;

imgFilename = sprintf( '%s', imgFilename); % remove blanks
imgFilename = [ imgFilename '.mat'];

try
   save( [ path imgFilename ], 'M');
catch
   err = 1;
end

cd( dir);