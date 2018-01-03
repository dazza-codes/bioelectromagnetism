% demonstration of laplacian interpolation

file = which('eeg_toolbox');
[eegpath, file, ext] = fileparts(file);
demoPath = fullfile(eegpath, 'lapint', '');
cd(demoPath)

if ~exist('pnt0012','var'),
  
  sphere_load;
  
  if isequal(exist('lapdemodata.mat'),2),
    load lapdemodata;
    % plot the dipole potentials at the surface
    figure
    subplot(2,2,1)
    patch('faces', tri0012,'vertices', pnt0012, 'facecolor', 'interp', 'Cdata', pot0012(:,1));
    title('12 vertices')
    subplot(2,2,2)
    patch('faces', tri0042,'vertices', pnt0042, 'facecolor', 'interp', 'Cdata', pot0042(:,1));
    title('42 vertices')
    subplot(2,2,3)
    patch('faces', tri0162,'vertices', pnt0162, 'facecolor', 'interp', 'Cdata', pot0162(:,1));
    title('162 vertices')
    subplot(2,2,4)
    patch('faces', tri0642,'vertices', pnt0642, 'facecolor', 'interp', 'Cdata', pot0642(:,1));
    title('642 vertices')
    fprintf('...loaded lapdemodata\n');
    return
  end
  
end


if ~exist('lap0012','var'),
  fprintf('...compute the laplacian matrix for the spherical triangulations\n');
	lap0012 = lapcal(pnt0012, tri0012);
	lap0042 = lapcal(pnt0042, tri0042);
	lap0162 = lapcal(pnt0162, tri0162);
	lap0642 = lapcal(pnt0642, tri0642);
	lap2562 = lapcal(pnt2562, tri2562);
end

if ~exist('lpot0012','var'),
  fprintf('...compute the laplacian of the potentials\n');
  lpot0012 = lap0012 * pot0012;
  lpot0042 = lap0042 * pot0042;
  lpot0162 = lap0162 * pot0162;
  lpot0642 = lap0642 * pot0642;
  lpot2562 = lap2562 * pot2562;
end

% determine interpolation matrices, interpolate all potentials
% and compute the residual variance for the interpolations

if ~exist('int0012_0042'),
	fprintf('...interpolate from 12 to 42/162/642/2562 vertices\n');
	int0012_0042 = lapint(lap0042, 1:12);	% interpolate from 12 to 42 vertices
	int0012_0162 = lapint(lap0162, 1:12);	% interpolate from 12 to 162 vertices
	int0012_0642 = lapint(lap0642, 1:12);	% interpolate from 12 to 642 vertices
	int0012_2562 = lapint(lap2562, 1:12);	% interpolate from 12 to 2562 vertices
end
if ~exist('rv0012_0042'),
  fprintf('...residual variance for 12 to 42/162/642/2562 vertices\n');
	rv0012_0042 = rv(pot0042, int0012_0042*pot0012);
	rv0012_0162 = rv(pot0162, int0012_0162*pot0012);
	rv0012_0642 = rv(pot0642, int0012_0642*pot0012);
	rv0012_2562 = rv(pot2562, int0012_2562*pot0012);
end

if ~exist('int0042_0162'),
	fprintf('...interpolate from 42 to 162/642/2562 vertices\n');
  int0042_0162 = lapint(lap0162, 1:42);	% interpolate from 42 to 162 vertices
	int0042_0642 = lapint(lap0642, 1:42);	% interpolate from 42 to 642 vertices
	int0042_2562 = lapint(lap2562, 1:42);	% interpolate from 42 to 2562 vertices
end
if ~exist('rv0042_0162'),
  fprintf('...residual variance for 42 to 162/642/2562 vertices\n');
	rv0042_0162 = rv(pot0162, int0042_0162*pot0042);
	rv0042_0642 = rv(pot0642, int0042_0642*pot0042);
	rv0042_2562 = rv(pot2562, int0042_2562*pot0042);
end

if ~exist('int0162_0642'),
  fprintf('...interpolate from 162 to 642/2562 vertices\n');
	int0162_0642 = lapint(lap0642, 1:162);		% interpolate from 162 to 642 vertices
	int0162_2562 = lapint(lap2562, 1:162);	% interpolate from 162 to 2562 vertices
end
if ~exist('rv0162_0642'),
  fprintf('...residual variance for 162 to 642/2562 vertices\n');
	rv0162_0642 = rv(pot0642, int0162_0642*pot0162);
	rv0162_2562 = rv(pot2562, int0162_2562*pot0162);
end

if ~exist('int0642_2562'),
  fprintf('...interpolate from 642 to 2562 vertices\n');
  int0642_2562 = lapint(lap2562, 1:642);
end
if ~exist('rv0642_2562'),
  fprintf('...residual variance for 642 to 2562 vertices\n');
  rv0642_2562 = rv(pot2562, int0642_2562*pot0642);
end

if ~isequal(exist('lapdemodata.mat'),2),
  save lapdemodata;
end

figure

subplot(2,2,1)
bar([rv0012_0042; rv0012_0162; rv0012_0642; rv0012_2562]')
title('interpolation from 12 vertices towards ...')
legend('42', '162', '642', '2562','Location','NorthWest')
axis([0 31 0 1]);

subplot(2,2,2)
bar([rv0042_0162; rv0042_0642; rv0042_2562]')
title('interpolation from 42 vertices towards ...')
legend('162', '642', '2562','Location','NorthWest')
axis([0 31 0 1]);

subplot(2,2,3)
bar([rv0162_0642; rv0162_2562]')
title('interpolation from 162 vertices towards ...')
legend('642', '2562','Location','NorthWest')
axis([0 31 0 1]);

subplot(2,2,4)
bar([rv0642_2562]')
title('interpolation from 642 vertices towards ...')
legend('2562','Location','NorthWest')
axis([0 31 0 1]);

return
