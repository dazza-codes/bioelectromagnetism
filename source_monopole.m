
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% monopole source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% a conductive medium
sigma = 0.2;

% with a point source current magnitude
Io = 10;

% located at the origin
IoR = [0,0,0];

% consider a spherical surface 
sphereR = 1;
FV = sphere_tri('ico',4,sphereR);

% By symmetry, the current flow is radially outward.
% At a distance r > 0 from the source, we have

% radial location, vector coordinates
R = FV.vertices;
R(:,1) = R(:,1) - IoR(1);
R(:,2) = R(:,2) - IoR(2);
R(:,3) = R(:,3) - IoR(3);
r = vector_magnitude(R,IoR);
r = repmat(r,1,3);

% extracellular currents at surface
Je = [ Io / (4 * pi) ] * [ R ./ (r.^3) ];

% electric potential at surface
V  = [ Io / (4 * pi * sigma) ] * [ 1 ./ r(:,1) ];

FV.Cdata = V;

figure;
H = patch('faces',FV.faces,'vertices',FV.vertices,...
  'FaceVertexCData',FV.Cdata,'facecolor','interp',...
  'edgecolor','none');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear FV

% a conductive medium
sigma = 0.2;

% with a point source current magnitude
Io = 10;

% located at the origin
IoR = [0,0,0];

FV.vertices = ones(10*10,3);
x = linspace(-10,10,10);
y = x;
n = 0;
for i = 1:length(x),
  for j = 1:length(y),
    n = n + 1;
    FV.vertices(n,:) = [ x(i), y(j), 10 ];
  end
end

FV.faces = delaunay(FV.vertices(:,1),FV.vertices(:,2));

% radial location, vector coordinates
R = FV.vertices;
R(:,1) = R(:,1) - IoR(1);
R(:,2) = R(:,2) - IoR(2);
R(:,3) = R(:,3) - IoR(3);
r = vector_magnitude(R,IoR);
r = repmat(r,1,3);

% extracellular currents at surface
Je = [ Io / (4 * pi) ] * [ R ./ (r.^3) ];

% electric potential at surface
V  = [ Io / (4 * pi * sigma) ] * [ 1 ./ r(:,1) ];

FV.Cdata = V;

figure;
H = patch('faces',FV.faces,'vertices',FV.vertices,...
  'FaceVertexCData',FV.Cdata,'facecolor','interp',...
  'edgecolor','none');
colormap(hot); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear FV

% a conductive medium
sigma = 0.2;

% with a point source current magnitude
Io = 10;

% located at the origin
IoR = [0,0,0];

dim = 20;

FV.vertices = ones(dim*dim*dim,3);
x = linspace(-10,10,dim);
y = x;
z = x;
n = 0;
for i = 1:length(x),
  for j = 1:length(y),
    for k = 1:length(z),
      n = n + 1;
      FV.vertices(n,:) = [ x(i), y(j), z(k) ];
    end
  end
end

FV.faces = delaunay3(FV.vertices(:,1),FV.vertices(:,2),FV.vertices(:,3));

% radial location, vector coordinates
R = FV.vertices;
R(:,1) = R(:,1) - IoR(1);
R(:,2) = R(:,2) - IoR(2);
R(:,3) = R(:,3) - IoR(3);
r = vector_magnitude(R,IoR);
r = repmat(r,1,3);

% extracellular currents at surface
Je = [ Io / (4 * pi) ] * [ R ./ (r.^3) ];

% electric potential at surface
V  = [ Io / (4 * pi * sigma) ] * [ 1 ./ r(:,1) ];


% create Analyze volume of potential
avw = avw_hdr_make;
avw.hdr.dime.dim(2:4) = [dim,dim,dim];

avw.img = zeros(dim,dim,dim);
n = 0;
for i = 1:length(x),
  for j = 1:length(y),
    for k = 1:length(z),
      n = n + 1;
      avw.img(i,j,k) = V(n);
    end
  end
end
avw_view(avw);
