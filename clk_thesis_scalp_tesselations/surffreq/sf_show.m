function fh = sf_show(data, varargin)
% SF_SHOW  	Show slices from 3D data
%   	    	fh = SF_SHOW(data [,colormap] [,fh] [,slices] [,val]);
%   	    	  COLORMAP = name of valid colormap
%   	    	  AH/FH = axis/figure handle to reuse (new figure by default)
%   	    	  SLICES = 1x3 array of indices for slices
%   	    	  VAL = threshold value to display

% Copyright (c) 1997 Clif Kussmaul. All rights reserved.

%%% NOTES

%%% THINGS TO DO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% action starts here

%%% check for valid arguments 
if nargin < 1
    help sf_show; return
end

%%% assign default values
cmap = 'gray';
ah = [];
fh = [];
val = 0;
x = size(data,1) / 2;
y = size(data,2) / 2;
z = size(data,3) / 2;

%%% interpret optional arguments
for j=1:nargin-1
    aj = varargin{j};
    if ischar(aj)
    	cmap = aj;
    elseif isempty(aj)
    elseif ishandle(aj) & strcmp(get(aj,'Type'),'axes')
    	ah = aj;
    elseif ishandle(aj) & strcmp(get(aj,'Type'),'figure')
    	fh = aj;
    elseif prod(size(aj)) == 1, 
    	val = aj;
    elseif all(size(aj) == [1 3])
    	x = aj(1); y = aj(2); z = aj(3);
    else
    	error(['Optional argument ', num2str(j), ' is invalid.']);
    end
    clear aj;
end

if any(size(data) == 1), data = squeeze(data); end
if ~isempty(ah)
    axes(ah);

elseif ~isempty(fh)
    set(0,'CurrentFigure',fh);
    
else
    fh = figure('Name', 	    	'sf_show', ...
	    	'PaperPositionMode',	'auto', ...
	    	'Position',		[0, 300, 500, 500]);
end
colormap(cmap);

%%% for 2D data, display image
if ndims(data) == 2
    set(0,'CurrentFigure',fh); subplot(1,1,1);
    if val
    	image(sf_find_3d(data, val, 'vol'), 'CDataMapping','scaled');
    else
    	image(data, 'CDataMapping','scaled');
    	sf_slice_fix(data);
    end

%%% for 3D data, display orthogonal slices
elseif ndims(data) == 3

    set(0,'CurrentFigure',fh); subplot(2,2,2);
    if val
    	image(sf_find_3d(squeeze(data(round(x),:,:)), val, 'vol'), ...
    	      'CDataMapping','scaled');
    else
    	image(squeeze(data(round(x),:,:)), 'CDataMapping','scaled');
        sf_slice_fix(data);
    end
    line([1 size(data,3)],[y y], 'LineWidth',1, 'Color','r');
    line([z z],[1 size(data,2)], 'LineWidth',1, 'Color','g');
    line([1 size(data,3) size(data,3)           1  1], ...
    	 [1           1  size(data,2) size(data,2) 1], ...
    	 'LineWidth',2, 'Color', 'b');
    title(['X=', num2str(x)]); xlabel('z'); ylabel('y');


    set(0,'CurrentFigure',fh); subplot(2,2,3);
    if val
    	image(sf_find_3d(squeeze(data(:,:,round(z))), val, 'vol'), ...
    	      'CDataMapping','scaled');
    else
    	image(squeeze(data(:,:,round(z))), 'CDataMapping','scaled');
        sf_slice_fix(data);
    end
    line([1 size(data,2)],[x x], 'LineWidth',1, 'Color','b');
    line([y y],[1 size(data,1)], 'LineWidth',1, 'Color','r');
    line([1 size(data,2) size(data,2)           1  1], ...
         [1           1  size(data,1) size(data,1) 1], ...
         'LineWidth',2, 'Color', 'g');
    title(['Z=', num2str(z)]); xlabel('y'); ylabel('x');


    set(0,'CurrentFigure',fh); subplot(2,2,4);
    if val
    	image(sf_find_3d(squeeze(data(:,round(y),:)), val, 'vol'), ...
    	      'CDataMapping','scaled');
    else
    	image(squeeze(data(:,round(y),:)), 'CDataMapping','scaled');
        sf_slice_fix(data);
    end
    line([1 size(data,3)],[x x], 'LineWidth',1, 'Color','b');
    line([z z],[1 size(data,1)], 'LineWidth',1, 'Color','g');
    line([1 size(data,3) size(data,3)           1  1], ...
    	 [1           1  size(data,1) size(data,1) 1], ...
    	 'LineWidth',2, 'Color', 'r');
    title(['Y=', num2str(y)]); xlabel('z'); ylabel('x');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% supporting functions start here

function sf_slice_fix(data)

if islogical(data)  	    caxis([0 1]);
elseif isa(data,'uint8')    caxis([0 255]);
else	    	    	    caxis([min(min(data)) max(max(data))]);
end
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% that's all folks

