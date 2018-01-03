function err = simpleIma2img( filename, fileList, path, params, hProgBar);
%
% Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
dir = pwd; cd( path);

% check filesize:
imaFid = fopen( char( fileList( 1,:)), 'r', 'ieee-be');
fseek( imaFid, 0, 'eof');
if ftell( imaFid) ~= 137216
   fprintf('\nSupport for images different to a 256*256 matrix not implemented !\n\n');
   fclose( imaFid);
   err = 1;
   return;
end
fclose( imaFid);

err = 0;
order = genSliceOrder( params.file.nSlices, params.user.verschachtelt);

imgFilename = [ filename '.img' ];
imgFid = fopen( imgFilename, 'w', 'ieee-le');


fprintf('processing: bitmaps\n');
for i = 1 : size(fileList,1); 
   imaFid = fopen( char( fileList( order(i),:)), 'r', 'ieee-be');  
   fseek( imaFid, 6144, 'bof');
   [img, n]= fread( imaFid, inf, 'uint16');
   if n ~= 65536
      fprintf('\nunable to read image completly. %2.4f %% readed => terminating !\n\n', n/65536);
      fclose('all');
      return;
   end
   m=fwrite( imgFid, img, 'uint16');
   if m ~= 65536
      fprintf('\nunable to write image => terminating !\n\n');
      fclose('all');
      return;
   end
   fclose( imaFid);
   if nargin == 5
      waitbar( i/size(fileList,1), hProgBar);
   end
end;
fclose(imgFid);

s = [params.file.matrix(1) params.file.matrix(2) size(fileList,1)];
v = [ (params.file.FoV(1) / params.file.matrix(1)) ...   
      (params.file.FoV(2) / params.file.matrix(2)) ...  
      (params.file.sliceThickness*(1+params.file.distFactor))];
spm_hwrite( filename, s, v, 1, 4, 0, [ 0 0 0], params.file.seqType);

% check carefully if img file exists before deleting the ima file
if params.user.delImas   
   try 
      fid = fopen( imgFilename, 'r');
      fopen( [imgFilename(1:end-4) '.hdr' ], 'r');
      fseek( fid, 0, 'eof');
      if ftell( fid) ~= (131072 * size(fileList,1)),
         err = 1;
      end
      fclose( 'all');
   catch
      err = 1;
   end
end

cd( dir);