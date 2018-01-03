function ImROIdemo()
%============================================================================
% demonstration of ImROI.m
%============================================================================
newline = sprintf('\n');
%============================================================================
% read or create an image
%============================================================================
I= imread('flowers.tif');
if isempty(I)
    I = peaks(128);
end;
currentfig = gcf;
if isempty(currentfig)
    currentfig = 0;
end
%============================================================================
% DEBUG test ROI statistics
%I = randn(128);
%============================================================================
% Run ImROI in different modes
%============================================================================
disp('Run ImROI manually drawing an ROI, without producing a file')
disp('roi = ImROI(I);')
disp(newline);
roi = ImROI(I);
%============================================================================
% Color the ROI as a patch
%============================================================================
disp('Using the output roi to create a patch with the average color');
disp('For 255 colordepth : meancolor = roi.mean/255; patch(roi.x,roi.y,meancolor;');
hold on

meancolor = roi.mean;
if length(meancolor) <3
    meancolor = meancolor(ones(1,3));
end;
meancolor = meancolor/255;  
patchhandle = patch(roi.x,roi.y,meancolor(1:3));
%============================================================================
% Wait to continue
%============================================================================
disp(newline);
disp('Hit return key to continue')
pause

%============================================================================
% Run ImROI with predefined x and y arrays, creating an output file
%============================================================================
disp(newline);
disp('Run ImROI with predefined x and y arrays')
[cx,cy] = size(I);
cx = floor(cx/2); cy = floor(cy/2);
r = min([cx,cy])/4;
phi = linspace(0,2*pi,72);
x = cx+r*cos(phi);y= cx+r*sin(phi);
disp('Create x and y arrays:')
disp('phi = linspace(0,2*pi,72)')
disp('x = cx+r*cos(phi);y= cx+r*sin(phi);')
disp('Command syntax: roi = ImROI(I,x,y, [])');
disp(newline)
disp('Enter a filename for output when prompted.');
disp('Enter a comment when prompted.');
roi = ImROI(I,x,y,[]);

%============================================================================
% Wait to exit
%============================================================================
disp('Hit return key to exit')
pause
%============================================================================
% Remove figure
%============================================================================
delete(currentfig+1:gcf);
return
%============================================================================
% end ImROIdemo
%============================================================================
