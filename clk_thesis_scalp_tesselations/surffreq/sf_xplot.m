function fh = sf_xplot(intr, extr, varargin)
% SF_XPLOT	Cross plots for each extrinsic coordinate
%		fh = SF_XPLOT(intr, extr [,axis] [,caxis] [,faces] [,fh] ...
%   	    	    	      [,interp] [,label] [,tmode] [,mode] )
%		plots INTRinsic and EXTRinsic points
%   	    	using 1x4  AXIS range for all axes
%   	    	using 1x2 CAXIS range for all color axes
%		using FACES or griddata as necessary
%   	    	  FH = figure handle to reuse (new figure by default)
%   	    	  INTERP = {'linear','cubic','nearest','invdist'}
%		  LABEL = string for figure title (appended)
%   	    	  TMODE = {'xyz','XYZ'}
%		  MODE = {'abs','mag','pha','magpha','real','imag','realimag'}

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES
% - FILL makes PATCH objects, PCOLOR makes SURFACE objects

%%% THINGS TO DO
% - add axis argument instead of global
% ? 3D results are wierd if one coordinate is constant
% ? avoid duplication of code for ex1, ex2
% ? color scale should not change within a movie...
% ? increase speed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments
if (nargin < 1) 
	help sf_xplot; return;
elseif (size(intr, 2) ~= 2)
    error(['Invalid number of columns (', num2str(size(intr, 2)), ') in INTR.']);
elseif all(size(extr, 2) ~= [1 3]) 
    error(['Invalid number of columns (', num2str(size(extr, 2)), ') in EXTR.']);
elseif (size(intr, 1) ~= size(extr, 1))
    error(['INTR and EXTR do not have equal number of rows (', ...
           num2str(size(intr, 1)), ' ~= ', num2str(size(extr, 1)), ').']);
end

%%% interpret optional arguments
  ar = [];
 car = [];
face = [];
  fh = [];
interp = 'invdist';
  mode = 'real';
 tmode = '';
    la = '';

for j = 1:nargin-2
    aj = varargin{j};
    if ischar(aj)
        switch aj
    	    case {'linear','cubic','nearest','invdist'}, 
    	    	interp = aj;
    	    case {'abs','mag','pha','magpha','real','imag','realimag'},
    	        mode = aj;
    	    case {'xyz','XYZ'},
    	    	tmode = aj;
	    otherwise,		
	    	la = [la, aj];
	end
    elseif isempty(aj)
    elseif ishandle(aj) & strcmp(get(aj,'Type'),'figure')
    	    	    	    	      fh = aj;
    elseif all(size(aj) == [1 2])    car = aj;
    elseif all(size(aj) == [1 4])     ar = aj;
    elseif size(aj,2) == 3	    face = aj;
    else
    	error(['Optional argument ', num2str(j), ' is invalid.']);
    end
    clear aj;
end

switch mode
    case 'abs',     	ex1 =   abs(extr);	ex2 = [];
    case 'mag',     	ex1 =   abs(extr);	ex2 = [];
    case 'pha',	    	ex1 = angle(extr);	ex2 = [];
    case 'magpha',	ex1 =   abs(extr);	ex2 = angle(extr);
    case 'real',	ex1 =  real(extr);	ex2 = [];
    case 'imag',	ex1 =  imag(extr);	ex2 = [];
    case 'realimag',	ex1 =  real(extr);	ex2 = imag(extr);
end

%%% divide into complex components

%%% get sizes and other key values
imin = min(min(intr));
imax = max(max(intr));
igap = 0.01 * (imax - imin);
imin = imin - igap;
imax = imax + igap;
emin = min(min(ex1));
emax = max(max(ex1));
if all(size(ex2))
    emin = min(min(min(ex2)),emin);
    emax = max(max(max(ex2)),emax);
end
numc = 2 * sqrt(length(intr));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for all plots 

if all(size(ex1)) & all(size(ex2))
    prow1 = 200;
    prow2 = 203;
    position = [ 0, 300, 300*size(extr,2), 500 ];
elseif all(size(ex1)) | all(size(ex2))
    prow1 = 100;
    prow2 = 100;
    position = [ 0, 300, 300*size(extr,2), 250 ];
end

prow1 = prow1 + 10 * size(extr, 2);
prow2 = prow2 + 10 * size(extr, 2);

if isempty(fh)
    fh = figure('Name', 			'sf_xplot', ...
	    	'PaperPositionMode',		'auto', ...
	    	'Position', 		    	position, ...
	    	'DefaultAxesDrawMode',	    	'fast', ...
	    	'DefaultLineLineWidth',	    	0.25, ...
	    	'DefaultLineMarkerSize',	2, ...
	    	'DefaultPatchEdgeColor',	'none', ...
	    	'DefaultPatchFaceColor',	'interp', ...
	    	'DefaultSurfaceEdgeColor',	'none', ...
	    	'DefaultSurfaceFaceColor',	'interp');
    sf_tool axes cmap misc shade zoom
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot using faces

if all(size(face))
    %%% use normals to remove backfaces
    ix = [intr, ones(size(intr,1),1)];
    fx = cross(ix(face(:,1),:) - ix(face(:,2),:), ...
    	       ix(face(:,1),:) - ix(face(:,3),:));
    fx = fx(:,3) < 0;

    for j = 1:size(ex1, 2)
    	set(0,'CurrentFigure',fh); subplot(prow1+j);
    	ax = gca;
    	patch('faces',	    	    	face(fx,:), ...
    	      'vertices',   	    	intr, ...
    	      'facevertexcdata',    	ex1(:, j), ...
      	      'facecolor',  	    	get(ax,'defaultsurfacefacecolor'), ...
              'edgecolor',  	    	get(ax,'defaultsurfaceedgecolor'));
    end
    for j = 1:size(ex2, 2)
    	set(0,'CurrentFigure',fh); subplot(prow2+j);
    	ax = gca;
    	patch('faces',	    	    	face(fx,:), ...
    	      'vertices',   	    	intr, ...
    	      'facevertexcdata',    	ex2(:, j), ...
      	      'facecolor',  	    	get(ax,'defaultsurfacefacecolor'), ...
              'edgecolor',  	    	get(ax,'defaultsurfaceedgecolor'));
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot using griddata

else
    [gu, gv] = meshgrid(linspace(imin, imax, numc));
    for j = 1:size(ex1,2)
	set(0,'CurrentFigure',fh); subplot(prow1+j); 
	pcolor(gu, gv, griddata(intr(:,1),intr(:,2), ex1(:,j), gu,gv, interp));
    end
    for j = 1:size(ex2,2)
	set(0,'CurrentFigure',fh); subplot(prow2+j); 
	pcolor(gu, gv, griddata(intr(:,1),intr(:,2), ex2(:,j), gu,gv, interp));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for all plots 

for j = 1:size(ex1,2)
    set(0,'CurrentFigure',fh); subplot(prow1+j);
    hold on
    plot(intr(:, 1), intr(:, 2), '.k'); 
    hold off
    axis('equal','square');
    if ~isempty(ar) axis(ar); end
    if ~isempty(tmode) title(tmode(j)); end
    if ~isempty(la) & (j == 2) title(la); end
    if isempty(car) caxis([emin,emax]), else caxis(car); end
    colorbar;
end
for j = 1:size(ex2,2)
    set(0,'CurrentFigure',fh); subplot(prow2+j);
    hold on
    plot(intr(:, 1), intr(:, 2), '.k'); 
    hold off
    axis('equal','square');
    if ~isempty(ar) axis(ar); end
    if ~isempty(tmode) title(tmode(j)); end
    if isempty(car) caxis([emin,emax]), else caxis(car); end
    colorbar;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

