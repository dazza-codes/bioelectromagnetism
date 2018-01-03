

% script to test vertex to voxel coordinate mapping


origin = avw_center(avw);

xo = origin.abs.mm(1);
yo = origin.abs.mm(2);
zo = origin.abs.mm(3);

% Find direction cosines for line from origin to nasion
x =   0 + xo;
y = 100 + yo;
z =   0 + zo;
r = sqrt( (x-xo).^2 + (y-yo).^2 + (z-zo).^2 );
l = (x-xo)/r; % cos alpha
m = (y-yo)/r; % cos beta
n = (z-zo)/r; % cos gamma

% find discrete points from origin to the vertex + 25%
POINTS = 250;
PERCENT = 0.25;
RMIN = 0;
rmax = r + (r * PERCENT);
radial_distances = linspace(RMIN,rmax,POINTS);

L = repmat(l,1,POINTS);
M = repmat(m,1,POINTS);
N = repmat(n,1,POINTS);

XI = (L .* radial_distances) + xo;
YI = (M .* radial_distances) + yo;
ZI = (N .* radial_distances) + zo;

% interpolate volume values at these points ( XI,YI are swapped here
% because the Analyze volume is a left handed coordinate system )
VI = interp3(vol,YI,XI,ZI,'*nearest');

% Find the intensity values in this array that approximate
% the scalp
