function [FV, Edges] = avw_shrinkwrap(avw,FV,interpVal,...
    fitval,fittol,fititer,fitchange,fitvattr)

% avw_shrinkwrap - Tesselate the surface of a 3D Analyze 7.5 avw struct
% 
% [FV, Edges] = avw_shrinkwrap(avw,FV,interpVal,...
%               fitval,fittol,fititer,fitchange,fitvattr)
% 
% avw       - an Analyze 7.5 data struct (see avw_read)
% FV        - input tesselation; if empty, sphere tesselation 
%             is created.  FV has fields FV.vertices, FV.faces
% interpVal - use radial interpolation (faster) or incremental
%             radial shrink (slower), 0|1 (default 1, faster)
% fitval    - image intensity to shrink wrap (default 20)
% fittol    - image intensity tolerance (default 5)
% fititer   - max number of iterations to fit (default 200)
% fitchange - least sig. change in intensity values 
%             between fit iterations (default 2)
% fitvattr  - vertex attraction (constraint), 0:1, smaller
%             values are less constraint; close to 0 for 
%             no constraint is useful when dealing with 
%             binary volumes, otherwise 0.4 (40%) seems OK
% 
% FV        - a struct with 2562 vertices and 5120 faces
% Edges     - a [2562,2562] matrix of edge connectivity for FV
% 
% This function has been developed to find the scalp surface 
% for MRI of the human head.  It is not a sophisticated, robust 
% algorithm!
% 
% It starts with a sphere tesselation (large radius) and moves
% each vertex point toward the center of the volume until it
% lies at or near the fitval.  The vertices are constrained to
% move only along the radial projection from the origin and they
% are also required to stay within a small distance of their
% neighbours.  The function is not optimised for speed, but 
% it should produce reasonable results.
% 
% Example of creating a scalp tesselation for SPM T1 MRI template:
% 
%   avw = avw_read('T1');
%   FV = avw_shrinkwrap(avw,[],0,0,[],intensity,5.0,50,0.5,0.4);
%   patch('vertices',FV.vertices,'faces',FV.faces,'facecolor',[.6 .6 .6]);
% 
% The output vertex coordinates are in meters with an origin at (0,0,0), 
% which lies at the center of the MRI volume.  This function uses the
% avw.hdr.dime.pixdim values to convert from voxel coordinates into
% meters.
% 
% See also: ISOSURFACE, SPHERE_TRI, MESH_REFINE_TRI4,
%           MESH_BEM_SHELLS_FUNC, MESH_BEM_SHELLS_SCRIPT
% 

% $Revision: 1.1 $ $Date: 2004/11/12 01:30:25 $

% Licence:  GNU GPL, no implied or express warranties
% History:  08/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - adapted from mesh_shrinkwrap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parse arguments

if ~exist('avw','var'), error('...no input volume\n');
elseif isempty(avw),    error('...empty input volume\n');
end

% Find approximate center of volume
center = avw_center(avw);
origin = center.abs.voxels;

version = '[$Revision: 1.1 $]';
fprintf('AVW_SHRINKWRAP [v%s]\n',version(12:16));  tic;

if ~exist('interpVal','var'), interpVal = 0;
elseif isempty(interpVal),    interpVal = 0;
end

if ~exist('fitval','var'),   fit.val = 20;
elseif isempty(fitval),      fit.val = 20;
else                         fit.val = fitval;
end

if ~exist('fittol','var'),   fit.tol = 5;
elseif isempty(fittol),      fit.tol = 5;
else                         fit.tol = fittol;
end

if fit.val <= fit.tol,
    error('...must use fit tolerance < fit value\n');
end

if ~exist('fititer','var'),  fit.iter = 200;
elseif isempty(fititer),     fit.iter = 200;
else                         fit.iter = fititer;
end

if ~exist('fitchange','var'),fit.change = 2;
elseif isempty(fitchange),   fit.change = 2;
else                         fit.change = fitchange;
end

if ~exist('fitvattr','var'), fit.vattr = 0.4;
elseif isempty(fitvattr),    fit.vattr = 0.4;
else                         fit.vattr = fitvattr;
end
if fit.vattr > 1,
    fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 1\n');
    fit.vattr = 1;
end
if fit.vattr < 0,
    fprintf('...fit vertattr (v) must be 0 <= v <= 1, setting v = 0.\n');
    fit.vattr = 0;
end


% get size of volume, in voxels
[xdim,ydim,zdim] = size(avw.img);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check whether to create a sphere tesselation
% or use an input tesselation as the start point

sphere = 0;
if ~exist('FV','var'),
    sphere = 1;
elseif ~isfield(FV,'vertices'),
    sphere = 1;
elseif ~isfield(FV,'faces'),
    sphere = 1;
elseif isempty(FV.vertices),
    sphere = 1;
elseif isempty(FV.faces),
    sphere = 1;
end
if sphere,
    % Create a sphere tesselation to encompass the volume
    diameter = max([xdim ydim zdim]);
    radius = diameter / 1.5;
    FV = sphere_tri('ico',4,radius); % 2562 vertices
    % Shift the center of the sphere to the center of the volume
    FV.vertices = FV.vertices + repmat(origin,size(FV.vertices,1),1);
else
    fprintf('...using input FV tesselation...\n');
end

% the 'edge' matrix is the connectivity of all vertices,
% used to find neighbours during movement of vertices,
% including smoothing the tesselation
FV.edge = mesh_edges(FV);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate scalp intensity
fit = estimate_scalp(FV,avw.img,origin,fit);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now begin recursion
fprintf('...fitting...\n');	tic;

i = 1;
Fminima = 0;
intensityAtRMean = [0 0];

while i <= fit.iter,
    
    if interpVal,
        % use radial interpolation method, moving directly
        % to the intensity value nearest correct intensity
        [FV, intensityAtR, R] = locate_val(FV,avw.img,origin,fit);
    else
        % use incremental method, moving along radial line
        % gradually until finding correct intensity
        [FV, intensityAtR, R] = shrink_wrap(FV,avw.img,origin,fit);
    end
    
    intensityAtRMean = [ intensityAtRMean(2), mean(intensityAtR) ] ;
    
    fprintf('...distance:  mean = %8.4f voxels, std = %8.4f voxels\n',mean(R),std(R));
    fprintf('...intensity: mean = %8.4f,        std = %8.4f\n',...
        mean(intensityAtR),std(intensityAtR));
    fprintf('...real iteration: %3d\n',i);
    
    % Is the mean distance reasonable?
    if mean(R) < 0.5,
        error('...mean distance < 0.5 voxel!\n');
    end
    
    % MDifVal is the mean of the absolute difference
    % between the vertex intensity and the fit intensity
    MDifVal = abs(intensityAtRMean(2) - fit.val);
    
    % Is the mean difference within the tolerance range?
    if MDifVal < fit.tol,
        fprintf('...mean intensity difference < tolerance (%5.2f +/- %5.2f)\n',...
            fit.val,fit.tol);
        break;
    else
        fprintf('...mean intensity difference > tolerance (%5.2f +/- %5.2f)\n',...
            fit.val,fit.tol);
    end
    
    % How much has the intensity values changed?
    if (i > 1) & intensityAtRMean(2) > 0,
        if intensityAtRMean(2) - intensityAtRMean(1) < fit.change,
            fprintf('...no significant intensity change (< %5.2f) in this iteration\n',...
                fit.change);
            Fminima = Fminima + 1;
            if Fminima >= 5,
                fprintf('...no significant intensity change in last 5 iterations\n');
                break;
            end
        else
            Fminima = 0;
        end
    end
    
    % Ensure that iterations begin when MDifVal is truly initialised
    if isnan(MDifVal),
        i = 1;
    else,
        i = i + 1;
    end
    
end

FV = mesh_smooth(FV,origin,fit.vattr);

% Remove large edges matrix from FV
Edges = FV.edge;
FV = struct('vertices',FV.vertices,'faces',FV.faces);

% Now center the output vertices at 0,0,0 by subtracting
% the volume centroid
FV.vertices(:,1) = FV.vertices(:,1) - origin(1);
FV.vertices(:,2) = FV.vertices(:,2) - origin(2);
FV.vertices(:,3) = FV.vertices(:,3) - origin(3);

% Scale the vertices by the pixdim values, after
% converting them from mm to meters
xpixdim = double(avw.hdr.dime.pixdim(2)) / 1000;
ypixdim = double(avw.hdr.dime.pixdim(3)) / 1000;
zpixdim = double(avw.hdr.dime.pixdim(4)) / 1000;

FV.vertices(:,1) = FV.vertices(:,1) .* xpixdim;
FV.vertices(:,2) = FV.vertices(:,2) .* ypixdim;
FV.vertices(:,3) = FV.vertices(:,3) .* zpixdim;


t=toc; fprintf('...done (%5.2f sec).\n\n',t);

return










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunctions...






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fit] = estimate_scalp(FV,vol,origin,fit),

xo = origin(1); yo = origin(2); zo = origin(3);

Nvert = size(FV.vertices,1);

% Estimate the scalp intensity, using 10% of vertices
N = round(0.05 * Nvert);
scalp_intensity = zeros(N,1);

fprintf('...estimating scalp intensity from %d vertices\n', N);
fprintf('...starting scalp intensity = %d\n', fit.val);

indices = round(rand(1,N) * Nvert);

i = 0;
for v = indices,
  
  x = FV.vertices(v,1);
  y = FV.vertices(v,2);
  z = FV.vertices(v,3);
  r = sqrt( (x-xo).^2 + (y-yo).^2 + (z-zo).^2 );
  l = (x-xo)/r; % cos alpha
  m = (y-yo)/r; % cos beta
  n = (z-zo)/r; % cos gamma
  
  % find discrete points from zero to the vertex
  POINTS = 250;
  radial_distances = linspace(0,r,POINTS);
  
  L = repmat(l,1,POINTS);
  M = repmat(m,1,POINTS);
  N = repmat(n,1,POINTS);
  
  XI = (L .* radial_distances) + xo;
  YI = (M .* radial_distances) + yo;
  ZI = (N .* radial_distances) + zo;
  
  % interpolate volume values at these points
  % ( not sure why have to swap XI,YI here )
  VI = interp3(vol,YI,XI,ZI,'*nearest');
  
  % find location in VI where intensity gradient is steep
  half_max = 0.5 * max(VI);
  index_maxima = find(VI > half_max);
  % use the largest index value to locate the maxima intensity that lie
  % furthest toward the outside of the head
  index_max = index_maxima(end);
  
  i = i + 1;
  scalp_intensity(i,1) = VI(index_max);
  
  plot(radial_distances,VI)
  
end

fit.val = mean(scalp_intensity);
fit.tol = std(scalp_intensity);

fprintf('...estimated scalp intensity = %g\n', fit.val);
fprintf('...estimated tolerance intensity = %g\n', fit.tol);

return










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, intensityAtR, R] = locate_val(FV,vol,origin,fit),

xo = origin(1); yo = origin(2); zo = origin(3);

Nvert = size(FV.vertices,1);
progress = round(.1 * Nvert);

% Initialise difference intensity & distance arrays
intensityAtR = zeros(Nvert,1);
R = intensityAtR;

% Find distance and direction cosines for line from 
% origin to all vertices
XV = FV.vertices(:,1);
YV = FV.vertices(:,2);
ZV = FV.vertices(:,3);
RV = sqrt( (XV-xo).^2 + (YV-yo).^2 + (ZV-zo).^2 );
LV = (XV-xo)./RV; % cos alpha
MV = (YV-yo)./RV; % cos beta
NV = (ZV-zo)./RV; % cos gamma

% Check for binary volume data, if empty, binary
binvol = find(vol > 1);

% Locate each vertex at a given fit value
tic
for v = 1:Nvert,
    
    if v > progress,
        fprintf('...interp3 processed %4d of %4d vertices',progress,Nvert);
        t = toc; fprintf(' (%5.2f sec)\n',t);
        progress = progress + progress;
    end
    
    % Find direction cosines for line from origin to vertex
    x = XV(v);
    y = YV(v);
    z = ZV(v);
    d = RV(v);
    l = LV(v); % cos alpha
    m = MV(v); % cos beta
    n = NV(v); % cos gamma
    
    % find discrete points between origin
    % and vertex + 20% of vertex distance
    POINTS = 250;
    
    Rarray = linspace(0,(d + .2 * d),POINTS);
    
    L = repmat(l,1,POINTS);
    M = repmat(m,1,POINTS);
    N = repmat(n,1,POINTS);
    
    XI = (L .* Rarray) + xo;
    YI = (M .* Rarray) + yo;
    ZI = (N .* Rarray) + zo;
    
    % interpolate volume values at these points
    % ( not sure why have to swap XI,YI here )
    VI = interp3(vol,YI,XI,ZI,'*linear');
    
    % do we have a binary volume (no values > 1)
    if isempty(binvol),
        maxindex = max(find(VI>0));
        if maxindex,
            R(v) = Rarray(maxindex);
        end
    else
        % find the finite values of VI
        index = max(find(VI(isfinite(VI))));
        if index,
            
            % Find nearest volume value to the required fit value
            nearest = max(find(VI >= fit.val));
            
            %[ nearest, value ] = NearestArrayPoint( VI, fit.val );
            
            % Check this nearest index against a differential
            % negative peak value
            %diffVI = diff(VI);
            %if max(VI) > 1,
            %    diffindex = find(diffVI < -20);
            %else
            % probably a binary volume
            %    diffindex = find(diffVI < 0);
            %end
            
            % now set d
            if nearest,
                R(v) = Rarray(nearest);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constrain relocation by fit.vattr, 
    % some % of distance from neighbours
    
    vi = find(FV.edge(v,:));  % the neighbours' indices
    X = FV.vertices(vi,1);    % the neighbours' vertices
    Y = FV.vertices(vi,2);
    Z = FV.vertices(vi,3);
    
    % Find neighbour distances
    RN = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
    % Find mean distance of neighbours
    neighbour_distances_mean = mean(RN);
    
    minattr = fit.vattr;
    maxattr = 1 + fit.vattr;
    
    if R(v) < (minattr * neighbour_distances_mean),
        R(v) = minattr * neighbour_distances_mean;
    end
    if R(v) > (maxattr * neighbour_distances_mean),
        R(v) = maxattr * neighbour_distances_mean;
    end
    if R(v) == 0, R(v) = neighbour_distances_mean; end
    
    % relocate vertex to new distance
    x = (l * R(v)) + xo;
    y = (m * R(v)) + yo;
    z = (n * R(v)) + zo;
    
    FV.vertices(v,:) = [ x y z ];
    
    % Find intensity value at this distance
    intensityAtR(v) = interp1(Rarray,VI,R(v),'linear');
    
end

if isempty(binvol),
    % check outliers and smooth twice for binary volumes
    FV = vertex_outliers(FV, R, origin);
    FV = mesh_smooth(FV,origin,fit.vattr);
end
FV = mesh_smooth(FV,origin,fit.vattr);

return









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV, intensityAtR, R] = shrink_wrap(FV,vol,origin,fit),

xo = origin(1); yo = origin(2); zo = origin(3);

Nvert = size(FV.vertices,1);
progress = round(.1 * Nvert);

% Initialise difference intensity & distance arrays
R = zeros(Nvert,1);
intensityAtR = R;

% Check for binary volume data, if empty, binary
binvol = find(vol > 1);

Imin = 0;
Imax = max(max(max(vol)));

% Manipulate each vertex
tic
for v = 1:Nvert,
  
  if v > progress,
    fprintf('...shrink-wrap processed %4d of %4d vertices',progress,Nvert);
    t = toc; fprintf(' (%5.2f sec)\n',t);
    progress = progress + progress;
  end
  
  % Find direction cosines for line from origin to vertex
  x = FV.vertices(v,1);
  y = FV.vertices(v,2);
  z = FV.vertices(v,3);
  r = sqrt( (x-xo).^2 + (y-yo).^2 + (z-zo).^2 );
  l = (x-xo)/r; % cos alpha
  m = (y-yo)/r; % cos beta
  n = (z-zo)/r; % cos gamma
  
  % interpolate volume values at this point ( x,y are swapped here
  % because the Analyze volume is a left handed coordinate system )
  intensity_old = interp3(vol,y,x,z,'*nearest');
  
  % move vertex closer to the origin, in a radial direction
  r_change = r * 0.05;
  r_new = r - r_change;
  
  % Calculate new vertex coordinates
  x = (l * r_new) + xo; % cos alpha
  y = (m * r_new) + yo; % cos beta
  z = (n * r_new) + zo; % cos gamma
  
  % interpolate volume values at this point ( x,y are swapped here
  % because the Analyze volume is a left handed coordinate system )
  intensity_new = interp3(vol,y,x,z,'*nearest');
  
  intensity_dif = intensity_new - intensity_old;
  
  if intensity_dif == 0,
    % relocate vertex to new distance
    FV.vertices(v,1) = x;
    FV.vertices(v,2) = y;
    FV.vertices(v,3) = z;
    intensityAtR(v,1) = intensity_new;
    R(v) = r_new;
  elseif (fit.val - intensity_new) < (fit.val - intensity_old),
    % relocate vertex to new distance
    FV.vertices(v,1) = x;
    FV.vertices(v,2) = y;
    FV.vertices(v,3) = z;
    intensityAtR(v,1) = intensity_new;
    R(v) = r_new;
  else
    intensityAtR(v,1) = intensity_old;
    R(v) = r;
  end
  
  FV = constrain_vertex(FV,v,origin);
  
end

FV = mesh_smooth(FV,origin,fit.vattr);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV] = constrain_vertex(FV,index,origin),

% This function adapts Smith, S. (2002), Fast robust automated brain
% extraction.  Human Brain Mapping, 17: 143-155.  It corresponds to update
% component 2: surface smoothness control.  It assumes that vertices are
% moved along a radial line from an origin, given by direction
% cosines, rather than calculating the surface normal vector.

xo = origin(1); yo = origin(2); zo = origin(3);

v = FV.vertices(index,:);
x = FV.vertices(index,1);
y = FV.vertices(index,2);
z = FV.vertices(index,3);

% Find radial distance of vertex from origin
r = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );

% Calculate unit vector
v_unit_vector = ( v - origin ) / r;

% Find direction cosines for line from center to vertex
l = (x-xo)/r; % cos alpha
m = (y-yo)/r; % cos beta
n = (z-zo)/r; % cos gamma

% Find neighbouring vertex coordinates
vi = find(FV.edge(index,:));  % the indices
neighbour_vertices = FV.vertices(vi,:);
X = neighbour_vertices(:,1);
Y = neighbour_vertices(:,2);
Z = neighbour_vertices(:,3);

% Find neighbour radial distances
r_neighbours = sqrt( (X-xo).^2 + (Y-yo).^2 + (Z-zo).^2 );
r_neighbours_mean = mean(r_neighbours);

% Find difference in radial distance between the vertex of interest and its
% neighbours; this value approximates the magnitude of sn in 
% Smith (2002, eq. 1 to 4)
r_diff = r - r_neighbours_mean;

% Find the vector sn, in the direction of the vertex of interest, given the
% difference in radial distance between vertex and mean of neighbours
sn = r_diff * v_unit_vector;

% Find distances between vertex and neighbours, using edge lengths.
% The mean value is l in Smith (2002, eq. 4)
edge_distance = FV.edge(index,vi);
edge_distance_mean = mean(edge_distance);

% Calculate radius of local curvature, solve Smith (2002, eq. 4)
if r_diff,
  radius_of_curvature = (edge_distance_mean ^ 2) / (2 * r_diff);
else
  radius_of_curvature = 10000;
end

% Define limits for radius of curvature
radius_min =  3.33; % mm
radius_max = 10.00; % mm

% Sigmoid function parameters,
% "where E and F control the scale and offset of the sigmoid"
E = mean([(1 / radius_min),  (1 / radius_max)]);
F = 6 * ( (1 / radius_min) - (1 / radius_max) );

Fsigmoid = (1 + tanh( F * (1 / radius_of_curvature - E))) / 2;

% multiply sigmoid function by sn
move_vector = Fsigmoid * sn;

FV.vertices(index,:) = v + move_vector;

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [val] = vol_val(vol,x,y,z),

% This function just ensures that xyz are
% actually within the volume before trying
% to get a volume value

val = nan; % assume zero value

x = round(x);
y = round(y);
z = round(z);

if x > 0 & y > 0 & z > 0,
    
    % get bounds of vol
    xdim = size(vol,1);
    ydim = size(vol,2);
    zdim = size(vol,3);
    
    if x <= xdim & y <= ydim & z <= zdim,
        % OK return volume value at xyz
        val = vol(x,y,z);
    end
end

return






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FV] = vertex_outliers(FV, R, origin),

xo = origin(1); yo = origin(2); zo = origin(3);

% Screen FV for outlying vertices, using
% mean +/- 2 * stdev of distance from origin
DistMean = mean(R);
DistStDev = std(R);
DistMax = DistMean + 2 * DistStDev;
DistMin = DistMean - 2 * DistStDev;

for v = 1:size(FV.vertices,1),
    
    if R(v) >= DistMax,
        R(v) = DistMean;
        relocate = 1;
    elseif R(v) <= DistMin,
        R(v) = DistMean;
        relocate = 1;
    else
        relocate = 0;
    end
    
    if relocate,
        x = FV.vertices(v,1);
        y = FV.vertices(v,2);
        z = FV.vertices(v,3);
        
        % Find direction cosines for line from center to vertex
        d = sqrt( (x-xo)^2 + (y-yo)^2 + (z-zo)^2 );
        l = (x-xo)/d; % cos alpha
        m = (y-yo)/d; % cos beta
        n = (z-zo)/d; % cos gamma
        
        % relocate vertex to this new distance
        x = (l * R(v)) + xo;
        y = (m * R(v)) + yo;
        z = (n * R(v)) + zo;
        
        FV.vertices(v,:) = [ x y z ];
    end
end
return
