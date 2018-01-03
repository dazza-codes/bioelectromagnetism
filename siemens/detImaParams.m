function params = detImaParams( filename);
%
% Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
% params = struct(
%                 name: string
%                 date: string
%                 time: string
%              seqType: string
%      acquisitionTime: double
%             normVect: [3x1 double]
%              colVect: [3x1 double]
%              rowVect: [3x1 double]
%          centerPoint: [3x1 double]
%    distFromIsoCenter: double
%                  FoV: [2x1 double]
%       sliceThickness: double
%           distFactor: double
%              repTime: double
%            scanIndex: double
%            angletype: [7x1 char]
%                angle: [4x1 char]
%               matrix: [2x1 douuble]
%              nSlices: double
fid = fopen( filename, 'r', 's');

% who and when
fseek( fid, 768, 'bof');
params.name = sscanf( char( fread( fid, 25, 'uchar')), '%s');
fseek( fid, 12, 'bof');
date = fread( fid, 3, 'uint32');
params.date = sprintf( '%d.%d.%d', date(3),date(2),date(1));
fseek( fid, 52, 'bof');
time = fread( fid, 3, 'uint32');
params.time = sprintf( '%d:%d:%d', time(1),time(2),time(3));

%scan stuff
fseek( fid, 3083, 'bof'); %sequenzeType
params.seqType = sscanf( char(fread( fid, 8, 'uchar')), '%s');

fseek( fid, 2048, 'bof'); % acquisition Time
params.acquisitionTime = fread( fid, 1, 'double');


%geometrical stuff
fseek( fid, 3792, 'bof');
params.normVect = fread( fid, 3, 'double');
fseek( fid, 3856, 'bof');
params.colVect = fread( fid, 3, 'double');
fseek( fid, 3832, 'bof');
params.rowVect = fread( fid, 3, 'double');
fseek( fid, 3768, 'bof');
params.centerPoint = fread( fid, 3, 'double');

fseek( fid, 3816, 'bof');
params.distFromIsoCenter = fread( fid, 1, 'double');

% sliceParams ...
fseek( fid, 3744, 'bof');
params.FoV = fread( fid, 2, 'double');

%fseek( fid, 5000, 'bof');
%params.pixelSize = fread( fid, 2, 'double');

fseek( fid, 1544, 'bof');
params.sliceThickness = fread( fid, 1, 'double');

switch params.seqType
case 'mpr'
   % here the number of slices seems to be at a different place ...   
   params.distFactor = 0;
   
otherwise
   %tested for mosaic-epi and t1-sag-screening seq.
   fseek( fid, 4136, 'bof');
   params.distFactor = fread( fid, 1, 'double');
end

   

fseek( fid, 1560, 'bof');
params.repTime = fread( fid, 1, 'double');

fseek( fid, 5726, 'bof');
params.scanIndex = str2num( sscanf( char(fread( fid, 3, 'uchar')), '%s'));

fseek( fid, 5814, 'bof');
params.angletype = char( fread( fid, 7, 'uchar'));
fseek( fid, 5821, 'bof');
params.angle = char( fread( fid, 4, 'uchar'));

fseek( fid, 5695, 'bof');
params.matrix(1) = str2num( sscanf( char(fread( fid, 3, 'uchar')), '%s'));
fseek( fid, 5700, 'bof');
params.matrix(2) = str2num( sscanf( char(fread( fid, 3, 'uchar')), '%s'));


%
switch params.seqType
case 'mpr'
   % here the number of slices seems to be at a different place ...   
   fseek( fid, 3984, 'bof');                       % if this fails try 3988 or 3992.
   params.nSlices = fread( fid, 1, 'uint32');
   
otherwise
   %tested for mosaic-epi and t1-sag-screening seq.
   fseek( fid, 4004, 'bof');
   params.nSlices = fread( fid, 1, 'uint32');
end

fclose( fid);