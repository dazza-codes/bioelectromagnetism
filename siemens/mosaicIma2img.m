function err = mosaicIma2img( fileList, path, params, hProgBar);
%
% Sebastian Thees 17.2.2001, email: s_thees@yahoo.com
%
% Dept. of Neurologie, Charite, Berlin, Germany
%
dir = pwd; cd( path);
%
% check filesize:
imaFid = fopen( char( fileList( 1,:)), 'r', 'ieee-be');

fseek( imaFid, 0, 'eof');

if params.file.nSlices < 16
   mosaicSize = 4;
elseif params.file.nSlices > 16 & params.file.nSlices < 64
   mosaicSize = 8;
else
   fprintf('\nSupport for images different to a 256*256 matrix not implemented !\n\n');
   fclose( imaFid);
   err = 1;
   return;
end
   
fclose( imaFid);

err = 0;
order = genSliceOrder( params.file.nSlices, params.user.verschachtelt);
%
fprintf('processing: bitmaps\n');

xSize = params.file.matrix(1);
ySize = params.file.matrix(2);
nX = mosaicSize;
nY = mosaicSize;
for n = 1 : mosaicSize^2, subimage{n}= uint16( zeros( xSize, ySize));   end

for i = 1 : size(fileList, 1);   
   imgFilename = sprintf( '%s', char(fileList(i,:)));   
   imgFilename = strrep( imgFilename, '.ima', '.img');
   imgFid = fopen( imgFilename, 'w', 'ieee-le');
   imaFid = fopen( char( fileList(i,:)), 'r', 'ieee-be');  
   fseek( imaFid, 6144, 'bof');   
   img = uint16( fread( imaFid, [mosaicSize*xSize mosaicSize*ySize], 'uint16'));
   if size( img, 1)*size(img,2) ~= mosaicSize*xSize*mosaicSize*ySize  % check for possible error
      fprintf('\nunable to read image completly. %2.4f %% readed => terminating !\n\n', ...
         size( img, 1)*size(img,2)/mosaicSize*xSize*mosaicSize*ySize);
      fclose('all');
      return;
   end
   %
   % begin: image extraction
   for j = 0 : nY - 1
      for k = 0 : nX - 1
         subimage{ 1+k+j*nX}(:,:) = img( 1 + k*xSize : (k+1)*xSize, 1 + j*xSize : (j+1)*xSize);
      end
   end
   % end: image extraction
   for n = 1 : params.file.nSlices
      m = fwrite( imgFid, subimage{ order(n)}, 'uint16'); % writing the subimages
      if m ~= xSize  * ySize % check for possible error
         err = 1;
      end
   end
   fclose( imaFid); fclose( imgFid);    
   s = [ xSize, ySize, params.file.nSlices];
   v = [ (params.file.FoV(1) / params.file.matrix(1)) ...   
         (params.file.FoV(2) / params.file.matrix(2)) ...  
         (params.file.sliceThickness*(1+params.file.distFactor))];
   spm_hwrite( imgFilename(1:end-4), s, v, 1, 4, 0, s/2, '');   
   
   % check carefully if img file exists before deleting the ima file
   if params.user.delImas   
      if err == 1
         fprintf('\nan error occured during img file creation, imas not deleted !\n\n');
         fclose( 'all');
         return;
      else
         try 
            fid = fopen( imgFilename, 'r', 'ieee-le');
            fopen( [imgFilename(1:end-4) '.hdr' ], 'r');
            fseek( fid, 0, 'eof');
            if ftell( fid) ~= (xSize*ySize * params.file.nSlices*2),
               err = 1;
            end
            fclose( 'all');
         catch
            err = 1;
         end
      end
   end
   if err == 1
      fprintf('\nan error occured during img file creation, terminating !\n\n');
      fclose( 'all');
      return;
   end;
   if nargin == 4
      waitbar( i/size(fileList,1), hProgBar);
   end
end;
cd( dir);

