function [ avw ] = ge_series2avw(examPath,seriesPath)

% ge_series2avw - converts a GE series to Analyze
% 
% avw = ge_series2avw(examPath,seriesPath)
% 
% Converts a series of GE slices into an Analyze 
% avw struct (see avw_read), which can be output
% as an Analyze .hdr/.img pair using avw_write.
% 
% examPath   - string path to an exam directory, 
%              which contains series directories 
%              below it
% seriesPath - the series to convert
%              (integer or string argument)
% 
% examPath is the name of the directory containing 
% the series subdirectories (e.g., series 1), which 
% contain the series image files (*.MR or I.*).  
% This function calls ge_series_read.
% 
% The function will attempt to reorient the GE 
% 3D volume into radiological orientation 
% (axial LAS, which is the default Analyze 
% orientation).  The resulting data should
% be SPM compatible when output with avw_write.
% 
% This function is in alpha development (as of 03/2003) 
% although a prior version has been tested with 
% Ax,Sag,Cor slices (with slice direction going both 
% ways). It was also tested for oblique axial, but 
% not on double obliques or anything more complicated.  
% The function does not provide information for an 
% SPM compatible .mat file.
% 
% see also ge_series_read, 
%          avw_view, avw_read, avw_write
%


% $Revision: 1.1 $ $Date: 2004/11/12 01:32:35 $

% Souheil J. Inati  <souheil.inati@nyu.edu> at 03/2003
% Dartmouth College, May 2000
% 
% Darren.Weber@flinders.edu.au, March 2003
% - Substantially redesigned file handling and function
%   call structures for integration with mri_toolbox at 
%   http://eeg.sf.net
% - Requested permission to distribute code under GPL licence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 2),
    doc ge_series2avw;
    error('...not enough input arguments.')
    return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the GE series header and image volume
[ge, lastfile] = ge_series_read(examPath, seriesPath);

% could try to use lastfile to create avw.fileprefix, but
% it is too variable to be reliable
% avw.fileprefix = lastfile

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert the GE series to an Analyze volume

% Generate the Analyze header
avw = ge_hdr2avw(ge);

version = '[$Revision: 1.1 $]';
fprintf('\nGE_SERIES2AVW [v%s]\n',version(12:16));  tic;

% Check if ADW scan (not sure this is useful, DLW 03/2003)
%if ge.hdr.image.user9 == 0, adwcount = 1;
%else,                    adwcount = ge.hdr.image.user9;
%end


% Reorient the GE data into radiological orientation during assignment
% of ge.img into avw.img (leave ge.img in original orientation)
avw = ge_reorient(ge, avw); % see below

% Write the Analyze files (not doing this)
%avw_write(avw,outName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the original write code, now replace with above (DLW)...

%outFile = strcat(outName,'.hdr');
%status = ge_writeSPMHeader(outFile,header);
%outFile = [outName sprintf('.img')];
%[fid,message] = fopen(outFile,'w');
%if (fid == -1),
%    fprintf('Cannot Open %s for writing.\n',outFile);
%    error(message);
%end
%fwrite(fid,reshape(imageVol,1,prod(size(imageVol))),'int16');
%status = fclose(fid);

t=toc; fprintf('...done (%5.2f sec).\n',t);

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [avw] = ge_reorient(ge, avw)

%ge_reorient - Assigns Analyze header dimensions and volume based on GE orientation
% 
% avw = ge_reorient(ge, avw)
% 
% reorients the GE 3D volume to be radiological
% orientation (axial LAS, which is SPM compatible)
% based on the GE acquition orientation
% 
% This has been tested with Ax,Sag,Cor with slices going
% both ways.  Also for Oblique axial. Don't count on double 
% obliques or anything really fancy.
% 


% Have looked over this orientation code carefully (DLW, 03/2003)
% The above comments are from a previous version, note it
% should be LAS, not RAS !!!!  I have found some inconsistencies
% in the code (03/2003) and fixed it as best I can for now.  Further
% testing with various volumes is required.


series_description = deblank(char(ge.hdr.series.se_desc)'); % unreliable!

fprintf('...checking GE series data orientation.\n');

% Determine the GE orientation
% orient is 1=axial, 2=sagittal, 3=coronal
% with opposite sign if backwards slice order
ras_start = char(ge.hdr.series.start_ras);
ras_end   = char(ge.hdr.series.end_ras);
if     strcmp(ras_start,'I') & strcmp(ras_end,'S'),
    fprintf('...slices are axial from inferior to superior.\n');
    orient =  1;
elseif strcmp(ras_start,'S') & strcmp(ras_end,'I'),
    fprintf('...slices are axial from superior to inferior.\n');
    orient = -1;
elseif strcmp(ras_start,'R') & strcmp(ras_end,'L'),
    fprintf('...slices are sagittal from right to left.\n');
    orient =  2;
elseif strcmp(ras_start,'L') & strcmp(ras_end,'R'),
    fprintf('...slices are sagittal from left to right.\n');
    orient = -2;
elseif strcmp(ras_start,'P') & strcmp(ras_end,'A'),
    fprintf('...slices are coronal from posterior to anterior.\n');
    orient =  3;
elseif strcmp(ras_start,'A') & strcmp(ras_end,'P'),
    fprintf('...slices are coronal from anterior to posterior.\n');
    orient = -3;
else,
    warning('GE orientation unknown!');
    orient = 0;
end

% Get the GE dimensions
nX = ge.hdr.image.imatrix_X;
nY = ge.hdr.image.imatrix_Y;
nZ = ge.hdr.image.slquant;     % slice quantity (Nslices)
pX = ge.hdr.image.pixsize_X;
pY = ge.hdr.image.pixsize_Y;
pZ = ge.hdr.image.slthick + ge.hdr.image.scanspacing;

[vX vY vZ] = size(ge.img);

% Reshape into axial radiological orientation (SPM compatible)
% The default Analyze orientation is +X left, +Y anterior, +Z superior (LAS)

switch orient,
    
case  0, % Unknown Orientation
    
    warning('avw.img = ge.img without reorientation!\n');
    avw.img = ge.img;
    
case {1, -1}, % Axial
    avw.hdr.dime.dim(2:4)    = [ nX nY nZ ];
    avw.hdr.dime.pixdim(2:4) = [ pX pY pZ ];
    
    avw.img = ge.img;
    if orient == 1, % Axial (I to S)
        % checked this (03/2003), not sure of L/R orient
        avw.img = flipdim(avw.img,2); % flip to P to A
    elseif orient == -1, % Axial (S to I)
        % have not checked this (03/2003)
        avw.img = flipdim(avw.img,2); % flip to P to A
        avw.img = flipdim(avw.img,3); % flip to I to S
    end
    
case {2, -2}, % Sagittal
    
    avw.hdr.dime.dim(2:4)    = [ nZ nX nY ];
    avw.hdr.dime.pixdim(2:4) = [ pZ pX pY ];
    
    avw.img = permute(ge.img,[3 1 2]);
    if orient == 2, % Sagittal (R to L)
        % have not checked this (03/2003)
        avw.img = flipdim(avw.img,2); % flip to P to A?
        avw.img = flipdim(avw.img,3); % flip to I to S?
    elseif orient == -2, % Sagittal (L to R)
        % checked this (03/2003)
        avw.img = flipdim(avw.img,1); % flip to R to L
        avw.img = flipdim(avw.img,2); % flip to P to A
        avw.img = flipdim(avw.img,3); % flip to I to S
    end
    
case {3, -3}, % Coronal
    
    avw.hdr.dime.dim(2:4)    = [ nX nZ nY ];
    avw.hdr.dime.pixdim(2:4) = [ pX pZ pY ];
    
    avw.img = permute(ge.img,[1 3 2]);
    if orient == 3, % Coronal (P to A)
        % have not checked this (03/2003), not sure of L/R orient
        avw.img = flipdim(avw.img,3); % flip to I to S?
    elseif orient == -3, % Coronal (A to P)
        % have not checked this (03/2003), not sure of L/R orient
        avw.img = flipdim(avw.img,2); % flip to P to A?
        avw.img = flipdim(avw.img,3); % flip to I to S?
    end
    
end

% Set the origin to the center of the volume (not sure this is valid, DLW)
%avw.hdr.dime.originator = [floor(avw.hdr.dime.dim(2)/2) ...
%                           floor(avw.hdr.dime.dim(3)/2) ...
%                           floor(avw.hdr.dime.dim(4)/2) 0 0];

avw.hdr.dime.glmax = max(max(max(avw.img)));
avw.hdr.hist.orient = '0';

return
