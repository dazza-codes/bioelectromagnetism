%%%%%%%%%%%%%%%%%%%%%%%%
%                      %
% Setup Analyze Header %
%                      %
%%%%%%%%%%%%%%%%%%%%%%%%
function [header, orient, im_offset, adwcount] = GE_createSPMHeader(fname)

% Read the info from the file
[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(fname);

% Initial ANALYZE Header Key
header.sizeOfHeader = uint32(348); 
header.dataType(1:10) = char(' ');
header.dbName(1:18) = char(' ');
header.extents = int32(16384);
header.sessionError = int16(0);
header.regular = char('r');
header.hkeyUnused = char(' ');
header.dim(1:8) = double(0);
%Header dim is actually INT16 but need to do arithmetic with it
header.voxUnits(1:4) = char(' ');
header.calUnits(1:8) = char(' ');
header.iUnused1 = int16(0);
header.datatype = int16(0);
header.bitpix = int16(0);
header.dimUnused = int16(0);
header.pixDim(1:8) = double(0);
header.voxOffset = double(0);
header.fUnused1 = double(0);
header.fUnused2 = double(0);
header.fUnused3 = double(0);
header.calMax = double(0);
header.calMin = double(0);
header.compressed = double(0);
header.verified = double(0);
header.glMax = int32(0);
header.glMin = int32(0);
header.description(1:80) = char(' ');
header.auxFile(1:24) = char(' ');
header.orient = int8(0);
header.originator(1:10) = char(' ');
header.generated(1:10) = char(' ');
header.scanNum(1:10) = char(' ');
header.patientID(1:10) = char(' ');
header.expDate(1:10) = char(' ');
header.expTime(1:10) = char(' ');
header.histUN0(1:3) = char(' ');
header.views = int32(0);
header.volsAdded = int32(0);
header.startField = int32(0);
header.fieldSkip = int32(0);
header.oMax = int32(0);
header.oMin = int32(0);
header.sMax = int32(0);
header.sMin = int32(0);

% Stuff in the file specific stuff
header.datatype = 4;
header.bitpix = 16;
header.orient = 0;
header.glMax = 65535;
header.dim(1) =  4;                % Dimensions (always 4?)
header.dim(2) =  im_hdr.imatrix_X; % X Voxels
header.dim(3) =  im_hdr.imatrix_Y; % Y Voxels
header.dim(4) =  im_hdr.slquant;   % Z Voxels
header.dim(5) =  1;                % Time Points
header.pixDim(2) = im_hdr.pixsize_X;                    % X Voxel Size
header.pixDim(3) = im_hdr.pixsize_Y;                    % Y Voxel Size
header.pixDim(4) = im_hdr.slthick + im_hdr.scanspacing; % Z Voxel Size


% Determine the orientation
if (strcmp(char(se_hdr.start_ras),'L'))  & (strcmp(char(se_hdr.end_ras),'R'))
     orient = 1;
elseif (strcmp(char(se_hdr.start_ras),'R'))  & (strcmp(char(se_hdr.end_ras),'L'))
     orient = -1;
elseif (strcmp(char(se_hdr.start_ras),'I'))  & (strcmp(char(se_hdr.end_ras),'S'))
     orient = 2;
elseif (strcmp(char(se_hdr.start_ras),'S'))  & (strcmp(char(se_hdr.end_ras),'I'))
     orient = -2;
elseif (strcmp(char(se_hdr.start_ras),'A'))  & (strcmp(char(se_hdr.end_ras),'P'))
     orient = 3;
elseif (strcmp(char(se_hdr.start_ras),'P'))  & (strcmp(char(se_hdr.end_ras),'A'))
     orient = -3;
else
     orient = 0;
end

% Check if ADW scan
if im_hdr.user9 == 0
  adwcount = 1;
else
  adwcount = im_hdr.user9;
end

return








