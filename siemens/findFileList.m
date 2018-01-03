function [fileList, n]= findFileList( filename, path, anatomOrIma)
%
% (c) Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
dir = pwd; cd( path);

ind1 = findstr( filename, '-');
ind2 = findstr( filename, '.ima');

switch anatomOrIma
case 1  % this is normal (anatom)
   prefix = filename( 1:ind1(2));
   fileNr = str2num( filename( ind1(2)+1 : ind2-1));

   file = filename;
   n    = 0;
   fid  = 0;
   while fid ~= -1
      fid = fopen( file, 'r');
      if fid ~= -1
         fclose( fid);
         n=n+1;
         file = [prefix num2str(fileNr-n) '.ima' ];
      end
   end
   m    = 0;
   fid  = 0;
   file = [prefix num2str(fileNr+1) '.ima' ];
   while fid ~= -1
      fid = fopen( file, 'r');
      if fid ~= -1
         fclose( fid);
         m=m+1;
         file = [prefix num2str(fileNr+m) '.ima' ];
      end
   end
   
   fileList = uint8( zeros( n+m-1, 1+size( [ prefix num2str(fileNr+m-1) '.ima' ], 2)));
   for i = 1 : n+m-1
      s = uint8( [ prefix num2str(fileNr-n+i) '.ima' ]);
      fileList(i,1:size(s,2)) = s;
   end
   n = n + m - 1;
   
case 2   % this means mosaic format (EPIs)
   prefix = filename( 1:ind1(1));
   fileNr(1) = str2num( filename( ind1(1)+1 : ind1(2)-1) );
   fileNr(2) = str2num( filename( ind1(2)+1 : ind2-1) );
   
   file = filename;
   n    = 0;
   fid  = 0;
   while fid ~= -1
      fid = fopen( file, 'r');
      if fid ~= -1
         fclose( fid);
         n=n+1;
         file = [prefix num2str(fileNr(1)-n) '-' num2str(fileNr(2)-n) '.ima' ];
      end
   end
   m    = 0;
   fid  = 0;
   file = [prefix num2str(fileNr(1)+1) '-' num2str(fileNr(2)+1) '.ima' ];
   while fid ~= -1
      fid = fopen( file, 'r');
      if fid ~= -1
         fclose( fid);
         m=m+1;
         file = [prefix num2str(fileNr(1)+m) '-' num2str(fileNr(2)+m) '.ima' ];
      end
   end
   if ~m, m=1; end
   fileList = uint8( zeros( n+m -1, 1+size( [ prefix num2str(fileNr(1)+m-1) '-' num2str(fileNr(2)+m-1) '.ima' ], 2)));
   for i = 1 : n+m-1
      s = uint8( [ prefix num2str(fileNr(1)-n+i) '-' num2str(fileNr(2)-n+i) '.ima' ]);
      fileList(i,1:size(s,2)) = s;
   end
   n = n + m - 1;
end

cd( dir);
