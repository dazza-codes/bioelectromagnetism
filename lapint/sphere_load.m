%Robert Oostenveld
%------------------------------------
%Department of Medical Physics      
%University Nijmegen, PO Box 9101
%The Netherlands, 6500 HB Nijmegen
%
% phone +31-24-3614230 or +31-24-3615213
% mail:roberto@mbfys.kun.nl
% http://www.mbfys.kun.nl/~roberto/
%------------------------------------

% I have computed the forward potential for a single sphere 
% using the BEM. The potential is computed on all vertices of
% each triangle, and I have created the sphere with different
% refinements.

% The original sphere12 is an icosahedron with 12 vertices, which 
% was stepwise refined to 42, 162, 642 and 2562 vertices. The 
% potential is computed on all vertices of the finest sphere.
% The file sphere42 is the first refinement of the icosahedron,
% and the first 12 vertices correspond to the vertices of 
% the file sphere12. This is the same for the other spheres, the 
% first vertices are always the ones corresponding to the coarser 
% meshes.

% The file z_axis2562.txt contains the potential. The format is 10*3 
% columns and 2562 rows, which represent the following:
% each row is the potential at one vertex
% column 1-3: dipole at z=0.0 with resp. x-, y- and z-orientation
% column 4-6: dipole at z=0.1 with resp. x-, y- and z-orientation
% column 7-9: dipole at z=0.2 with resp. x-, y- and z-orientation
% etc.
% I think that colums 1:3:30 (in Matlab notation) are the interesting ones 
% for interpolation, they are a tangential source which moves from the
% center to the surface.

% Since I used a regular refinement, you can take a selection of the 
% potential matrix (the transfer matrix) to get the potential for 
% each of the spheres. The first 12 electrodes correspond with sphere12,
% the first 42 correspond with sphere42, etc.

% Therefore, interpolating the potential of the first 12 electrodes (rows)
% of the potential (which correspond with the first 12 vertices of
% sphere42) onto all 42 vertices should result in a potential similar to 
% the first 42 rows of the potential matrix. You can repeat this for all
% combinations of triangulated spheres.

if isequal(exist('spheredata.mat'),2),
    load spheredata;
    return
end


% read the triangulated surfaces, each is a sphere of radius 1.0
[pnt0012, tri0012] = read_tri('sphere0012.tri');
[pnt0042, tri0042] = read_tri('sphere0042.tri');
[pnt0162, tri0162] = read_tri('sphere0162.tri');
[pnt0642, tri0642] = read_tri('sphere0642.tri');
[pnt2562, tri2562] = read_tri('sphere2562.tri');

% read the potential for a dipole moving along the z-axis
load z_axis2562.txt

% select the corresponding channels from the complete transfer matrix
% for each triangulation
pot0012 = z_axis2562(1:12,:);
pot0042 = z_axis2562(1:42,:);
pot0162 = z_axis2562(1:162,:);
pot0642 = z_axis2562(1:642,:);
pot2562 = z_axis2562(1:2562,:);

%figure; triplot(pnt0162, tri0162, pot0162(:, (10-1)*3+1),'contour'); rotate3D on


save spheredata

return


% make some 3D plots showing the data
for i=1:10,
    col = (i-1)*3;
    % plot the potential for a dipole oriented in x-direction
    subplot(10,3,col+1);
    triplot(pnt0162, tri0162, pot0162(:, col+1),'surface');
    % plot the potential for a dipole oriented in y-direction
    subplot(10,3,col+2);
    triplot(pnt0162, tri0162, pot0162(:, col+2),'surface');
    % plot the potential for a dipole oriented in z-direction
    subplot(10,3,col+3);
    triplot(pnt0162, tri0162, pot0162(:, col+3),'surface');
end
