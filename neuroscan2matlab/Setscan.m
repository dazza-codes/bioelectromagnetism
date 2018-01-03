%
%SETSCAN Workspace for NeuroScan and Macro Utilities.
%	SETSCAN Set up some variables needed with utilities. Should be run at 
%	the beginning of session to avoid errors, use startup.
%
%DIAGNOSTICS
%	Used global variables are R0...R9 OLD_R, S0...S9, UIHELP_FIG,
%	UIHELP_TXT, UIHELP_FILES, PATHS, SHOWWAIT_TXT, BRK, DIALOGX, DIALOGY,
%	BORDERX, BORDERY and BORDERB. Needs DIRS to be predefined global.
%
%EXAMPLES
%	Setscan

%Mention source when using or modifying these Shareware tools

%JVIR, 8-Apr-1999 Added some help text.
%JVIR,Jussi.Virkkala@occuphealth.fi

%JVIR, 4-Feb-1999 By default showing menubars.
%JVIR, 3-Feb-1999 Increased BORDERY
%JVIR, 2-Feb-1999 Modified for PCWIN Matlab 5.2.

%	J.Virkkala 19-Aug-94
%	J.Virkkala 30-Dec-94 
%	J.Virkkala  6-Mar-95 Part of ScanUtil
%	J.Virkkala 18-May-95 Together setutil and setscan.

disp(' ');
disp('setscan  : NeuroScan and Macro Utilities 2.1:')

disp('         : -clear figures.')
c=get(0,'children');
while ~isempty(c);close(c(1));c=get(0,'children');end
close			% close figures

disp('         : -global variables.')
global R0 R1 R2 R3 R4 R5 R6 R7 R8 R9 OLD_R 
global S0 S1 S2 S3 S4 S5 S6 S7 S8 S9
global UIHELP_FIG UIHELP_TXT UIHELP_FILE
global SHOWWAIT_TEXT BRK 

global DIALOGX DIALOGY BORDERX BORDERY BORDERB
disp('         : -dialog box position.')
pos=get(0,'screensize');
DIALOGX=pos(3)/2;DIALOGY=20;
BORDERX=8;
BORDERY=40; %20 before menus, 12
BORDERB=50;


BRK=0;

disp('         : -default figure properties.')
set(0,'units','pixels');
%JVIR, set(0,'defaultfiguremenubar','none')	% large figures
set(0,'defaultfigurepapertype','a4')
set(0,'defaultfigurepaperunits','centim')
%set(0,'defaultfigurepaperposition',[1 2 18 27])
set(0,'defaultfigureunits','pixels')	
set(0,'defaultaxesfontsize',10)
set(0,'defaulttextfontsize',10)

disp('         : -path for data and calibration files.');
disp(' ');
disp('examples : qeega1');
global SCANPATH
%JVIR
%SCANPATH=findpath('scanutil');
SCANPATH=[cd DIRS]; 
path(path,[SCANPATH 'Headers']);
path(path,[SCANPATH 'Files']);
path(path,[SCANPATH 'Help']);
path(path,[SCANPATH 'Functions']);
disp(' ');

%END OF SETSCAN