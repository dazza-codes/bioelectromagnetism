

% source dipole script

eegPath = eeg_toolbox_path;
cd(eegPath)

LMAX = 250;

if exist('source_dipole.mat') == 2,
  load source_dipole
else
  Cdata = cell(LMAX,1);
  for L = 1:LMAX,
    FV = source_dipole_sphere([],[],L);
    Cdata(L) = FV.Cdata;
  end
  FV.Cdata = Cdata; clear Cdata;
  save source_dipole FV;
end

Nvertices = 6;
vertices = round(rand(Nvertices,1) * 1000);
vertex_potentials = zeros(Nvertices,LMAX);

for Lmax = 1:LMAX,
  
  vertex_potentials(:,Lmax) = FV.Cdata{Lmax}(vertices);
  
  if LMAX < 15,
    figure;
    patch('vertices',FV.vertices,'faces',FV.faces,'FaceVertexCData',FV.Cdata{Lmax},'Facecolor','interp');
    colorbar; view([-40,35])
  end
end

LPLOT = 1:20;

COLS = 3;
ROWS = ceil(Nvertices / COLS);

figure
PLOT = 0;
for v = 1:Nvertices,
  PLOT = PLOT + 1;
  subplot(ROWS,COLS,PLOT);
  plot(LPLOT,vertex_potentials(v,LPLOT))
  hold on
  scatter(LPLOT,vertex_potentials(v,LPLOT))
end
