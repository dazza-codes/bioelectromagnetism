
% plots_topography - script to plot ERP topo maps

clear

comp = 'wm';
data = 'link14hz';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set time for topo maps, in msec
xmin = 50;
xmax = 805;
xstep = 5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the data files
switch comp,
case 'sa',
    exp = 'oac'; con = 'ouc'; dif = 'oac-ouc';
case 'wm',
    exp = 'tac'; con = 'oac'; dif = 'tac-oac';
case 'ea',
    exp = 'oat'; con = 'oac'; dif = 'oat-oac';
case 'dt',
    exp = 'tud'; con = 'tuc'; dif = 'tud-tuc';
end
plots_loaddata


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output directories
if mkdir('topomaps'),
    cd('topomaps');
    if mkdir(comp),
        cd(comp);
    else
        msg = sprintf('Cannot create %s directory',comp);
        error(msg);
    end
else
    error('Cannot create ''topomaps'' directory');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup the eeg toolbox parameters

p = eeg_toolbox_defaults;
p.volt.sampleHz = 400;
p.volt.sampleMsec = 2.5;
p.volt.channels = 124;
p.volt.epochStart = -200;
p.volt.epochEnd = 1500;
p.volt.sweeps = 1;
p.volt.interpZero = 1;

% define time point 0
p.volt.sampleTime = 0;
p.volt.samplePoint = 81;
p.timeMethod = 1;
p.endTime = 81;

% Setup to plot the scalp mesh topography
p.elec.plot = 0;
p.elec.plotSurf = 0;
p.mesh.plotSurf = 1;

% Setup the contours
p.contour.plot3D = 1; % with contours

% Open the default mesh and electrodes
p = mesh_open(p);
p = elec_open(p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create topography maps for the exp/con conditions

expvolt = eval(CONT.expvolt);
convolt = eval(CONT.convolt);
CONTexpvolt = expvolt;
CONTconvolt = convolt;

expvolt = eval(PTSD.expvolt);
convolt = eval(PTSD.convolt);
PTSDexpvolt = expvolt;
PTSDconvolt = convolt;

% -- Determine common scale

CONT.erpmax = max([CONTexpvolt CONTconvolt]);
CONT.erpmin = min([CONTexpvolt CONTconvolt]);
CONT.ymax = ceil( max(CONT.erpmax));
CONT.ymin = floor(min(CONT.erpmin));

PTSD.erpmax = max([PTSDexpvolt PTSDconvolt]);
PTSD.erpmin = min([PTSDexpvolt PTSDconvolt]);
PTSD.ymax = ceil( max(PTSD.erpmax));
PTSD.ymin = floor(min(PTSD.erpmin));

ymin = min(CONT.ymin,PTSD.ymin);
ymax = max(CONT.ymax,PTSD.ymax);

ymax = max(abs(ymin),abs(ymax));
ymin = -1 * ymax;






% -- Load the first voltage data file

p.volt.path = datapath;
p.volt.file = CONT.expfile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf






% -- Load the second voltage data file

p.volt.path = datapath;
p.volt.file = CONT.confile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf









% -- Load the third voltage data file

p.volt.path = datapath;
p.volt.file = PTSD.expfile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf









% -- Load the fourth voltage data file

p.volt.path = datapath;
p.volt.file = PTSD.confile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create topography maps for the difference data

CONTdifvolt = eval(CONT.difvolt);

PTSDdifvolt = eval(PTSD.difvolt);

% -- Determine common scale

CONT.difmax = max(CONTdifvolt);
CONT.difmin = min(CONTdifvolt);
CONT.ymax = ceil( max(CONT.difmax));
CONT.ymin = floor(min(CONT.difmin));

PTSD.difmax = max(PTSDdifvolt);
PTSD.difmin = min(PTSDdifvolt);
PTSD.ymax = ceil( max(PTSD.difmax));
PTSD.ymin = floor(min(PTSD.difmin));

ymax = ceil( max([CONT.ymax PTSD.ymax]));
ymin = floor(min([CONT.ymin PTSD.ymin]));

ymax = max(abs(ymin),abs(ymax));
ymin = -1 * ymax;




% -- Load the first dif voltage data file

p.volt.path = datapath;
p.volt.file = CONT.diffile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf







% -- Load the second dif voltage data file

p.volt.path = datapath;
p.volt.file = PTSD.diffile;
p = eeg_open(p);

% -- Plot the topography

p = eeg_contours_engine(p);

% -- Setup the Animation parameters
H = get(gcf,'userdata');

% no mouse rotation
set(H.Info,'value',0);

% specify ymin/ymax
set(H.Yset,'value',1);
set(H.Ymin,'string',sprintf('%7.2f',ymin));
set(H.Ymin,'value',ymin);
set(H.Ymax,'string',sprintf('%7.2f',ymax));
set(H.Ymax,'value',ymax);

% specify x start/finish/step
set(H.Start,'value',xmin);
set(H.Start,'string',sprintf('%7.2f',xmin));
set(H.Finish,'value',xmax);
set(H.Finish,'string',sprintf('%7.2f',xmax));
set(H.Step,'value',xstep);
set(H.Step,'string',sprintf('%7.2f',xstep));

% specify save graphics
set(H.Save,'value',1);

set(H.gui,'userdata',H);

% -- For each view (Front,Back,Left,Right)
for view = 3:6,
    
    set(H.View,'value',view); mouse_rotate('view',H);
    set(H.gui,'userdata',H);
    gui_topo_animate('animate');
    
end

close gcf




return
