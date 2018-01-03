function [XYZ] = freesurfer_crs2xyz(CRS)


T = [ -1  0  0  128;
       0  0 -1 -128;
       0 -1  0  128 ];

size(CRS')

XYZ = [T * CRS]';

return
