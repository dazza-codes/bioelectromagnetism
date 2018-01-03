% script openFlatWindow
%
% Sets up the FLAT data structure, opens and initializes the
% Flat window.
%
% djh and baw, 7/98
%
% Modifications:
%
% djh, 4/99
% - Eliminate overlayClip sliders
% - Added mapWin sliders to show overlay only for pixels with parameter
%   map values that are in the appropriate range.
% wap & rfd, 9/99
% - Allow FLAT.subdir to be changed/set elsewhere
% bw 12/29/00
% - scan slider instead of buttons; support for types started

global mrSESSION

cd(mrSESSION.homeDir);

% Close Flat window if it exists
if exist('FLAT','var')
   if isfield(FLAT,'subdir')
      subdir = FLAT.subdir;	% FLAT might get obliterated, so save this
   end
   if ~isempty(FLAT) & isfield(FLAT,'ui')
		if isfield(FLAT.ui, 'windowHandle'), close(FLAT.ui.windowHandle); end
   end
end

disp('Initializing Flat view')

global FLAT

FLAT.name='FLAT';
FLAT.type='FLAT';
if exist('subdir','var')
	FLAT.subdir=subdir;
else
	FLAT.subdir='Flat';
end

% Refresh function, gets called by refreshScreen
FLAT.refreshFn = 'refreshView';

%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize data slots %
%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize slot for anat
FLAT.anat = [];

% Initialize slots for co, amp, and ph
FLAT.co = [];
FLAT.amp = [];
FLAT.ph = [];
FLAT.map = [];

% Initialize slots for tSeries
FLAT.tSeries = [];
FLAT.tSeriesScan = NaN;
FLAT.tSeriesSlice = NaN;

% Initialize ROIs
FLAT.ROIs = [];
FLAT.selectedROI = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute/load FLAT.coords and FLAT.grayCoords %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FLAT = getFlatCoords(FLAT);

%%%%%%%%%%%%%%%%%%%%%%%%
% Compute FLAT.ui.mask %
%%%%%%%%%%%%%%%%%%%%%%%%

FLAT = makeFlatMask(FLAT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make displayModes and color maps %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize displayModes
FLAT=resetDisplayModes(FLAT);
FLAT.ui.displayMode='anat';

%%%%%%%%%%%%%%%%%%%
% Open the window %
%%%%%%%%%%%%%%%%%%%

% Figure number for flat view window
FLAT.ui.figNum=figure('MenuBar','none');

% Handle for flat view window
FLAT.ui.windowHandle = gcf;

% Handle for main axis of flat view
FLAT.ui.mainAxisHandle = gca;

% Adjust position and size of main axes so there's room for
% buttons, colorbar, and sliders
set(FLAT.ui.mainAxisHandle,'position',[0.025 0.1 0.8 0.7]);

% Set window title
set(gcf,'Name',sprintf('Flat: %s',mrSESSION.homeDir));

% Set minColormap property so there's potentially room for 128
% colors 
set(FLAT.ui.windowHandle,'minColormap',128)

% Sharing of colors seems like it might be OK, but I'm turning it
% off just to be sure (djh, 1/26/98).
set(gcf,'sharecolors','off');

% Set closeRequestFcn so we can clean up when the window is closed
set(gcf,'CloseRequestFcn','closeFlatWindow');

%%%%%%%%%%%%%
% Add Menus %
%%%%%%%%%%%%%

disp('Attaching flat menus')

FLAT=filesMenu(FLAT);
FLAT=windowMenu(FLAT);
FLAT=analysisFlatMenu(FLAT);
FLAT=viewMenu(FLAT); 
FLAT=roiMenu(FLAT);
FLAT=plotMenu(FLAT); 
FLAT=colorMenu(FLAT);
FLAT=xformFlatMenu(FLAT);

%%%%%%%%%%%%%%%%%
% Add Color Bar %
%%%%%%%%%%%%%%%%%

% Make color bar and initialize it to 'off'
FLAT.ui.colorbarHandle=makeColorBar(FLAT);
setColorBar(FLAT,'off');
FLAT.ui.cbarRange = [];

%%%%%%%%%%%%%%%
% Add Buttons %
%%%%%%%%%%%%%%%

disp('Attaching buttons')

% Make buttons for choosing hemisphere
FLAT=makeHemisphereButtons(FLAT);
setCurSlice(FLAT,1);

%%%%%%%%%%%%%%%
% Add sliders %
%%%%%%%%%%%%%%%

disp('Attaching sliders')

% Scan number slider
w = 0.12; h = 0.04; l = 0.45; b = 0.95;
FLAT = makeSlider(FLAT ,'scan',[1,mrSESSION.nScans],[l b w h]);
FLAT = initScanSlider(FLAT ,1);

% correlation threshold:
FLAT = makeSlider(FLAT,'cothresh',[0,1],[.85,.85,.15,.05]);
setCothresh(FLAT,0);

% phase window:
FLAT = makeSlider(FLAT,'phWinMin',[0,2*pi],[.85,.75,.15,.05]);
FLAT = makeSlider(FLAT,'phWinMax',[0,2*pi],[.85,.65,.15,.05]);
setPhWindow(FLAT,[0 2*pi]);

% parameter map window: 
FLAT = makeSlider(FLAT,'mapWinMin',[0,1],[.85,.55,.15,.05]);
FLAT = makeSlider(FLAT,'mapWinMax',[0,1],[.85,.45,.15,.05]);
setMapWindow(FLAT,[0 1]);

% anatClip: determines clipping of the anatomy base-image
%           values to fill the range of available grayscales.
FLAT = makeSlider(FLAT,'anatMin',[0,1],[.85,.2,.15,.05]);
FLAT = makeSlider(FLAT,'anatMax',[0,1],[.85,.1,.15,.05]);
setAnatClip(FLAT,[0 .5]);

%%%%%%%%%%%%%%%%%%%
% Add popup menus %
%%%%%%%%%%%%%%%%%%%

disp('Attaching popup menus')
FLAT = makeTypePopup(FLAT,'Type','original',[.85,.95,.15,.05]);

FLAT = makePopup(FLAT,'ROI','none',[.85,.325,.15,.05]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize display parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize image field
FLAT.ui.image = [];

% Show all ROIs
FLAT.ui.showROIs = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%
% Load user preferences %
%%%%%%%%%%%%%%%%%%%%%%%%%

FLAT = loadPrefs(FLAT);

%%%%%%%%%%%%%%%%%%
% Refresh screen %
%%%%%%%%%%%%%%%%%%

FLAT=refreshScreen(FLAT);

disp('Done initializing Flat view')

