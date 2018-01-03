% NeuroScan and Macro Utilities toolbox.
% Version 2.1 1999-Apr-8 Small modifications.
% Version 2.0 10-Feb-1999 For PCWIN Matlab 5.2.
% Copyright (c) 1994-1999 Jussi Virkkala. 
% For PCWIN Matlab 5.2. 
% Functions to be executed: startup.
%
% NeuroScan data translation.
%   loadcnt    - Load continuous NeuroScan data.
%   scanhead   - General part of the NeuroScan header.
%   scanelec   - Electrode specific part of the NeuroScan header.
%   scaneven   - Eventtable of NeuroScan file.
%
% Display data.
%   axes1020   - Return axes handles for each electrode position.
%   convtext   - Convert plottext output into separate text objects.
%   getband    - Get bands from axes.
%   grid1020   - Interpolate values between electrodes.
%   makeblt    - Prepare current figures for gray screen capture.
%   plotband   - Plots vertical bands with multiple values.
%   plotdata   - Plot data with labels into current axes.   
%   plottext   - Plot text with numerical values into current axes.
%   pos1020    - Plot scalp electrode positions in XY or XYZ.
%   scrollda   - To scroll data, with uiscroll.
%   uiscroll   - Add slide control to x-axes.
%
% Analyze data.
%   fftspect   - Calculate power, amplitude spectrum using FFT.
%
% UserInterface.
%   claa       - Clears figure's axes inside normalized box.
%   creafig    - Create figure at wanted location.
%   creafigs   - Create figures spacely around screen.
%   crearefe   - Create LaTex reference.
%   creahtml   - Create HTML reference, new for 2.0
%   errorr     - Error response.
%   o2s        - Return handles as strings.
%   figtext    - Create figure text. 
%   showwait   - Show calculation time on figures name.
%   startup    - Startup for Utilities.
%   mutexc     - Create mutually exclusive uicontrols.
%   selelec    - Select EEG electrodes from scalp.
%   seltext    - Select strings.
%   setscan    - Workspace for NeuroScan Utilities functions.
%   setobj     - Set object properties.
%   startup    - Startup for utilities.
%   uihelp     - Read ascii string and shows it in help figure.
%   uicput     - Create uicontrols.
%   uicsetup   - Setup basic uicontrols at the top of figure.
%   uicval     - Return string which return objects property.
%
% Developing macros.
%   editmacr   - Create figure to load, edit and save macros.
%   evalstr    - Evaluate string as m file.
%   loadedit   - Load ascii string for editmacr and scanedit.
%   makem      - Create m files name0..9.m from R0..R9 macro.
%   maker      - Create macro R0..R9 from m files name0..9.m.
%
% General extensions to matlab functions.
%   deblanks   - Strip blanks in the beginning and end.
%   filename   - Return name from fid
%   findpath   - Return given path from MATLAB directory search path.
%   fpos       - Reposition fid.
%   freeze     - Freeze/Unfreeze figure by toggling visible.
%   loadstr    - Load ascii file.
%   putstr     - Put rows of string into another string.
%   rmrow      - Remove indicated rows from matrice.
%   savestr    - Save string as ascii into file. 
%   scrimage   - Display handle as image.
%   sortstr    - Sort string in ascending order.
%   strfind    - Extension to findstr, accept matrices.
%   tiestr     - Tie rows of strings into one matrice.
%   untiestr   - Untie strings from matrice.
%
% Example files.
%   demorgb    - Demonstrate RGB colors, using uicontrols.
%   scantut0   - Create tutorial text files.                            
%   uicview    - Add uicontrols to fig to control view in 3D. 
%   qeega1     - Quantitative EEG analysis, start this.
%   qeega2     - Quantitative EEG analysis, electrode plotting.
%   qeega3     - Quantitative EEG analysis, calculation.
%
% Auxiliary files*
%   Header     - Directory, saved calibration values.
%   Files      - Directory, data files.
%   Help       - Directory, help and documentation files.
%   Results    - Directory, result files.
%

%Mention source when using or modifying these tools
%JVIR, Jussi.Virkkala@occuphealth.fi
%JVIR,  9-Feb-1999 Changed directories.
%JVIR,  2-Feb-1999 Modified for PCWIN Matlab 5.2.

% 	J.Virkkala  3-Apr-95 Part of ScanUtil.
%	J.Virkkala 18-May-95 Added some functions to make toolbox complete.
%	J.Virkkala  7-Jun-95 Version for B.Sc.

%END OF CONTENTS                        		