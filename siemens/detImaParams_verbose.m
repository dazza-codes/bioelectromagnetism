function params = detImaParams_verbose();

global d;

[filename path] = uigetfile( '*.ima');

fid = fopen( [path filename], 'r', 's');

% who and when
fseek( fid, 768, 'bof');
params.name = char( fread( fid, 25, 'uchar'));
fseek( fid, 12, 'bof');
date = fread( fid, 3, 'uint32');
params.date = sprintf( '%d.%d.%d', date(3),date(2),date(1));
fseek( fid, 52, 'bof');
time = fread( fid, 3, 'uint32');
params.time = sprintf( '%d:%d:%d', time(1),time(2),time(3));
fprintf( '\n\n\n\n\nImage of    >>%s<<   recorded %s at %s.\n\n', params.name, params.date, params.time );


%sequenceType
fseek( fid, 3083, 'bof');
params.seqType = sscanf( char(fread( fid, 8, 'uchar')), '%s');
fprintf( 'sequence:  >>%s<<\n', params.seqType);

fseek( fid, 2048, 'bof'); % acquisition Time
params.acquisitionTime = fread( fid, 1, 'double');
fprintf( 'acquisitionTime [s]:  %4.8f\n\n', params.acquisitionTime);



%geometrical stuff
fseek( fid, 3792, 'bof');
params.normVect = fread( fid, 3, 'double');
fprintf( 'normVect (Ex):    Z=%4.8f     Y=%4.8f    X=%4.8f\n', params.normVect(3), params.normVect(2), params.normVect(1));
fseek( fid, 3856, 'bof');
params.colVect = fread( fid, 3, 'double');
fprintf( 'colVect (Ey):     Z=%4.8f     Y=%4.8f    X=%4.8f\n', params.colVect(3), params.colVect(2), params.colVect(1));
fseek( fid, 3832, 'bof');
params.rowVect = fread( fid, 3, 'double');
fprintf( 'rowVect (Ez):     Z=%4.8f     Y=%4.8f    X=%4.8f\n', params.rowVect(3), params.rowVect(2), params.rowVect(1));
fseek( fid, 3768, 'bof');
params.centerPoint = fread( fid, 3, 'double');
fprintf( 'centerPoint:      Z=%4.8f     Y=%4.8f    X=%4.8f\n\n', params.centerPoint(3), params.centerPoint(2), params.centerPoint(1));

fseek( fid, 3816, 'bof');
params.distFromIsoCenter = fread( fid, 1, 'double');
fprintf( 'distance from isocenter:  %4.2f\n', params.distFromIsoCenter);
fprintf( 'abs(normalVector)  = %4.8f\n', sqrt( params.normVect'*params.normVect ) );
fprintf( 'abs(centerPoint)   = %4.8f\n\n', sqrt( params.centerPoint' * params.centerPoint ) );

% sliceParams ...
fseek( fid, 3744, 'bof');
params.FoV = fread( fid, 2, 'double');
fprintf( 'FOV:        row=%4.8f   col=%4.8f\n', params.FoV(1), params.FoV(2));

%fseek( fid, 5000, 'bof');
%params.pixelSize = fread( fid, 2, 'double');
%fprintf( 'PixelSize:  row=%4.8f    col=%4.8f\n', params.pixelSize(1), params.pixelSize(2) );


fseek( fid, 1544, 'bof');
params.sliceThickness = fread( fid, 1, 'double');
fprintf( 'sliceThickness:  %4.8f\n', params.sliceThickness);

fseek( fid, 4136, 'bof');
params.distFactor = fread( fid, 1, 'double');
fprintf( 'distFactor:  %4.8f\n', params.distFactor);

fseek( fid, 1560, 'bof');
params.repTime = fread( fid, 1, 'double');
fprintf( 'repititionTime:  %4.8f\n', params.repTime);

fseek( fid, 5726, 'bof');
params.scanIndex = str2num( sscanf( char(fread( fid, 3, 'uchar')), '%s'));
fprintf( 'scanIndex:  %d\n\n', params.scanIndex);

fseek( fid, 5814, 'bof');
params.angletype = char( fread( fid, 7, 'uchar'));
fseek( fid, 5821, 'bof');
params.angle = char( fread( fid, 4, 'uchar'));
fprintf( 'angle: %s, size: %s.\n\n', params.angletype, params.angle);
%
switch params.seqType
case 'mpr'
   % here the number of slices seems to be at a different place ...   
   fseek( fid, 3984, 'bof');                       % if this fails try 3988 or 3992.
   params.nSlices = fread( fid, 1, 'uint32');
   fprintf( 'nSlices:  %4d\n', params.nSlices);
   
otherwise
   %tested for mosaic-epi and t1-sag-screening seq.
   fseek( fid, 4004, 'bof');
   params.nSlices = fread( fid, 1, 'uint32');
   fprintf( 'nSlices:  %4d\n', params.nSlices);  
end
