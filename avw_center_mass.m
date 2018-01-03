function center = avw_center_mass(avw)

% avw_center_mass - find center of mass for an image volume
%
% [center] = avw_center_mass(avw)
%
% avw  - an Analyze 7.5 data struct, see avw_read
%
% center is the Cartesian coordinates for
% the volume center of mass.  It is a struct with
% coordinates in both voxels and mm (1x3).
% 
% see also avw_center
%

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no implied or express warranties
% History:  09/2004, Darren.Weber_at_radiology.ucsf.edu
%                    Tania Boniske
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

version = '[$Revision: 1.1 $]';
fprintf('\nAVW_CENTER_MASS [v%s]\n',version(12:16));  tic;

% Extract info from avw.hdr
xdim = double(avw.hdr.dime.dim(2));
ydim = double(avw.hdr.dime.dim(3));
zdim = double(avw.hdr.dime.dim(4));

xpix = double(avw.hdr.dime.pixdim(2));
ypix = double(avw.hdr.dime.pixdim(3));
zpix = double(avw.hdr.dime.pixdim(4));

% For 1 dimensional space, the center of gravity is:
%
% x_pos = ( sum(i = 1:n) ( mi * xi ) ) / sum(i = 1:n) mi;
% y_pos = ( sum(i = 1:n) ( mi * yi ) ) / sum(i = 1:n) mi;
% z_pos = ( sum(i = 1:n) ( mi * zi ) ) / sum(i = 1:n) mi;
%
% If mi = 1, then the above is simply the mean(xi) = 1/n * sum(xi).
%
% note how the sum of mi is in both the numberator/denominator, so the
% units of the weights cancel out, leaving only the spatial units of xi.
%
% For 2 dimensions, the mi becomes 2D, so mi = sum(1:i,1:j) m(i,j);
% For 3 dimensions, the mi becomes 3D, so mi = sum(1:i,1:j,1:k) m(i,j,k);
%
% mx(i) is intended to caluculate the scalar yz mass at each x(i).
% xpos is then the division of the sum of the product of these scalar masses by the
% position by the sum of the scalar masses. sum(mass*position)/sum(mass).
% if the mass = 1 then the COG is the average position.

xi = 1:xdim;
yj = 1:ydim;
zk = 1:zdim;

mx = zeros(xdim,1);
my = zeros(ydim,1);
mz = zeros(zdim,1);

for i = xi,
    mx(i,1) = sum(sum(avw.img(i,:,:)));
end
xpos = (xi*mx)/sum(mx);

for j = yj,
    my(j,1) = sum(sum(avw.img(:,j,:)));
end
ypos = (yj*my)/sum(my);

for k = zk,
    mz(k,1) = sum(sum(avw.img(:,:,k)));
end
zpos = (zk*mz)/sum(mz);

% Find center voxel of volume
center.voxels = [xpos, ypos, zpos];

% Find mm coordinates of that voxel
center.mm = center.voxels .* [xpix, ypix, zpix];

t =toc; fprintf('...done (%6.2f sec)\n\n',t);

return
