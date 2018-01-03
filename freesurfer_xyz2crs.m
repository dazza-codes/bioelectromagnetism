function [XYZ] = freesurfer_xyz2crs(CRS)


T = [ -1  0  0  128;
       0  0 -1 -128;
       0 -1  0  128;
       0  0  0    1 ];


one = ones(size(CRS,1),1);

CRS = [CRS';one']';

