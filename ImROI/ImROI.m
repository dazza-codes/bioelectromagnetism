function roi = ImROI(im, x, y , outfilename, plotopt)
%============================================================================
% roi = ImROI(im, x, y , outfilename, plotopt)
% filename        : ImROI.m
% directory       :  
% description     : Calculates the area and mean pixel value and SD
%			        in a closed patch defined by the vectors x,y 
%			        (first and last points have to be connected)
%
% parameters in   : image data 
% and optional    : x,y         : x, y vector of points defining the patch
%                   outfilename : string with name for output filename 
%                                 Enter empty argument to choose a file
%                                 Data will be appended!
%		           : plotopt    : switches plotting of the image 
%                                 on (1) or off (0) default = 1 (on)
%
% parameters out  : roi structure
% 		            with roi x and y,roi statistics and timestamp
% Example 
% I = imread('flowers.tif');
% roi = ImROI(I,[],[],[])
% run ImROIdemo.m
%============================================================================
linecolor = 'r';
%============================================================================
% Display the image and get image size
%============================================================================
fighandle = figure;
imhandle = imagesc(im);
axis image; axis off;
axishandle = gca;

[nrows,ncols] = size(im);
%============================================================================
% take ROI x and y arrays from input arguments
% forece x,y lenghts to be equal and truncate values at image boundaries
%============================================================================
if (nargin < 5) 
    %============================================================================
    % plotopt = 0 suppresses the ROI plot
    %============================================================================
    plotopt = 1;
end;
if (nargin < 4) 
  outfilename = '';
  fid = 0; 
elseif length(outfilename) > 0 
  fid = fopen(outfilename, 'a');
elseif isempty(outfilename)
  [outfilename, outpathname] = uiputfile('*', 'Select an output file');
  fid = fopen([outpathname outfilename ], 'a');
  if fid < 0
    errmsg = sprintf('Warning: could not open file %s for append',...
                        outfilename);
    disp(errmsg)
  end
end;

if (nargin >= 3) 
    if ~isempty(x) & ~isempty(y) 
        nx = length(x);
        ny = length(y);
        x= x(1:min(nx,ny));
        y= y(1:min(nx,ny));  
        x( x < 0.5 ) = 0.5;
        x( x > ncols+0.5 ) = ncols+0.5;
        y( y < 0.5 ) = 0.5;
        y( y > nrows+0.5 ) = nrows+0.5;
        hold on; 
        linehandle = plot(x,y, linecolor);
    end;
else
   %============================================================================
   % Get the ROI interactively
   %============================================================================
   [x , y, linehandle] = getpoints(axishandle);
end;

%============================================================================
%Calculate ROI area
%============================================================================
n = length(x);
diffx = [diff(x) (x(1) - x(n)) ];
diffy = [diff(y) (y(1) - y(n)) ];
avector = y .* diffx + diffx .* diffy ./2;

%============================================================================
% Copy area, vectors and linehandle to roi stucture
%============================================================================
roi.apix = abs(sum(avector));
roi.x = x;
roi.y = y; 
roi.linehandle = linehandle;

%============================================================================
% Change the pointer to something that is familiar to Microsoft users...
%============================================================================
oldpointershape = get(fighandle, 'Pointer');
set(fighandle, 'Pointer', 'watch');

%============================================================================
%Calculate the ROI area in square mm
%============================================================================
XData = get(imhandle, 'XData'); 
YData = get(imhandle, 'YData');
pixarea = (diff(XData) +1) * (diff(YData) + 1);
   
%============================================================================
% Create the smallest rectangular grid around the ROI
%============================================================================
xmingrid = max( XData(1), floor(min(x))  );
xmaxgrid = min( XData(2),  ceil(max(x))  );
ymingrid = max( YData(1), floor(min(y))  );
ymaxgrid = min( YData(2),  ceil(max(y))  );  
xgrid = xmingrid:xmaxgrid;
ygrid = ymingrid:ymaxgrid;
[X, Y] = meshgrid(xgrid, ygrid);
cdata = get(imhandle, 'CData');
smallcdata = double(cdata(ygrid,xgrid,:));
[m,n,ncolors] = size(smallcdata);
%============================================================================
% Analyze only the points in the polygon
%============================================================================
k_inside= inpolygon(X,Y, x,y);
Xin = X(k_inside); 
Yin = Y(k_inside);
clear X Y
%============================================================================
% Determine the center of the polygon
%============================================================================
roi.center =  [mean(Xin(:)),mean(Yin(:))];
clear Xin Yin
%============================================================================
% Calculate the mean, SD, etc... and as fields add to roi structure
% do this for each color
%============================================================================
for i=1:ncolors
    roicidata = smallcdata(:,:,i);
    roi.mean(i) =   mean(roicidata(k_inside));
    roi.std(i)  =    std(roicidata(k_inside));
    roi.min(i)  =    min(roicidata(k_inside));
    roi.max(i)  =    max(roicidata(k_inside)); 
    roi.median(i)= median(roicidata(k_inside));
end;
%============================================================================
% Add the date and time for future reference (in files)
%============================================================================
roi.timestamp = datestr(now);
%============================================================================
% Remove line if plotopt = 0
%============================================================================
if ~plotopt
    delete(linehandle)
end;
%============================================================================
% Reset pointer shape
%============================================================================
set(fighandle, 'Pointer', oldpointershape);
      
%============================================================================
% Output to screen and to file
%============================================================================
  
  
if fid > 0
    comment = input('Comment : ', 's');
    disp(['Writing to file ', outfilename])
    fprintf(fid, '%-20s, %-20s\n', 'timestamp', roi.timestamp );    
    fprintf(fid, '%-20s, %6.2f\n', 'pix area', roi.apix);  
    fprintf(fid, '%-20s, %6.2f,%6.2f\n', 'roicenter', roi.center); 
    fprintf(fid, '%-20s, %-40s\n','comment', comment); 
    
    fprintf(fid, '%-20s, ', 'mean');  
    fprintf(fid, '%6.2f,', roi.mean);  
    fprintf(fid, '\n');  
    fprintf(fid, '%-20s, ', 'std');  
    fprintf(fid, '%6.2f,', roi.std);  
    fprintf(fid, '\n');  
    fprintf(fid, '%-20s, ', 'min');  
    fprintf(fid, '%6.2f,', roi.min);  
    fprintf(fid, '\n');  
    fprintf(fid, '%-20s, ', 'max');  
    fprintf(fid, '%6.2f,', roi.max);  
    fprintf(fid, '\n');  
    fprintf(fid, '%-20s, ', 'median');  
    fprintf(fid, '%6.2f,', roi.median);  
    fprintf(fid, '\n');  
  
    fprintf(fid, '\n');  
    fclose(fid)
else
    comment = [];
end; 

disp(sprintf('%-20s, %-20s', 'timestamp', roi.timestamp) );
disp(sprintf('%-20s, %6.2f, %6.2f', 'roicenter [x,y]', roi.center));  
disp(sprintf('%-20s, %6.2f', 'pix area', roi.apix) ); 
if ~isempty(comment)
    disp(sprintf('%-20s, %-40s',  'comment', comment) ); 
end;
disp(sprintf('%-20s, %6.2f,%6.2f, %6.2f', 'mean', roi.mean) );  
disp(sprintf('%-20s, %6.2f,%6.2f, %6.2f', 'std',  roi.std) );  
disp(sprintf('%-20s, %6.2f,%6.2f, %6.2f', 'min', roi.min) );  
disp(sprintf('%-20s, %6.2f,%6.2f, %6.2f', 'max', roi.max) );   
disp(sprintf('%-20s, %6.2f,%6.2f, %6.2f', 'median', roi.median) );     
%============================================================================
% LOCAL FUNCTION GETPOINTS
%============================================================================
function [xs,ys, linehandle] = getpoints(axishandle)

%============================================================================
% Find parent figure for the argument axishandle
%============================================================================
axes(axishandle);
figure(get(axishandle, 'Parent'));
%===========================================================================
% Change pointer shape
%===========================================================================
oldpointershape = get(gcf,'Pointer');

ptrc =  ones(16)+1;
ptrc( 1, :) = 1; 
ptrc(16, :) = 1; 
ptrc(: , 1) = 1; 
ptrc(: ,16) = 1; 
ptrc(1:4,8:9) = 1;
ptrc(8:9,1:4) = 1;
ptrc(13:16, 8:9 ) = 1;
ptrc( 8:9 ,13:16) = 1;
ptrc(5:12,5:12) = NaN;
set(gcf,'Pointer', 		'custom',...
 	'PointerShapeCData', 	ptrc,...
 	 'PointerShapeHotSpot',	[8 8]);

%===========================================================================
% Prepare for interactive collection of ROI boundary points
%===========================================================================
hold on
pointhandles = [];
xpts = [];
ypts = [];
splinehandle= [];
n = 0;
but = 1;
BUTN = 0;
KEYB = 1;
done =0;
%===========================================================================
% Loop until right hand mouse button or keayboard is pressed
%===========================================================================
while ~done;  
  %===========================================================================
  % Analyze each buttonpressed event
  %===========================================================================
  keyb_or_butn = waitforbuttonpress;
  if keyb_or_butn == BUTN;
    currpt = get(axishandle, 'CurrentPoint');
    seltype = get(gcf,'SelectionType');
    switch seltype 
    case 'normal',
      but = 1;
    case 'alt',
      but = 2;
    otherwise,
      but = 2;
    end;          
  elseif keyb_or_butn == KEYB
    but = 2;
  end; 
  %===========================================================================
  % Get coordinates of the last buttonpressed event
  %===========================================================================
  xi = currpt(2,1);
  yi = currpt(2,2);
  %===========================================================================
  % Start a spline throught the points or 
  % update the line through the points with a new spline
  %===========================================================================
  if but ==1
        if ~isempty(splinehandle)
           delete(splinehandle);
        end;
    	pointhandles(n+1) = plot(xi,yi,'ro');
	n = n+1;
	xpts(n,1) = xi;
	ypts(n,1) = yi;
	%===========================================================================
	% Draw a spline line through the points
    %===========================================================================
	if n > 1
	  t = 1:n;
	  ts = 1: 0.1 : n;
	  xs = spline(t, xpts, ts);
	  ys = spline(t, ypts, ts);
	  splinehandle = plot(xs,ys,'r-');
	end;
  elseif but > 1
      %===========================================================================
	  % Exit for right hand mouse button or keyboard input
      %===========================================================================
      done = 1;
  end;
end;

%===========================================================================
% Add first point to the end of the vector for spline 
%===========================================================================
xpts(n+1,1) = xpts(1,1);
ypts(n+1,1) = ypts(1,1);

%===========================================================================
% (re)draw the final spline 
%===========================================================================
if ~ isempty(splinehandle)
    delete(splinehandle);
end;    

t = 1:n+1;
ts = 1: 0.25 : n+1;
xs = spline(t, xpts, ts);
ys = spline(t, ypts, ts);

linehandle = plot(xs,ys,'r-');
drawnow;
%===========================================================================
% Delete the point markers 
%===========================================================================
if ~isempty(pointhandles)
    delete(pointhandles)
end;
%===========================================================================
% Reset pointershape 
%===========================================================================
set(gcf,'Pointer',oldpointershape);

%===========================================================================
% END OF LOCAL FUNCTION GETPOINTS 
%===========================================================================
