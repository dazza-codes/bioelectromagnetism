function BRIK2MGH(fname,xdim,ydim,zdim,timepoints,type,foutstem)
% BRIK2MGH function to read 3D+time AFNI BRIK
% of any type and dimensionality and convert
% to MGH-style bshort of bfloat files and 
% accompanying MGH header files
%
% Example: this writes out 16 .bshort and .hdr files
% BRIK2MGH('/data9/mike/AS/afni/ASt3avvr+orig.BRIK',64,64,16,100,'short','test_')
%
% Written December 11, 1998 by Timothy M. Ellmore
% Laboratory of Brain and Cognition, NIMH

disp(' ')
disp('**********************************')
disp(' ')
disp('BRIK2MGH: reading raw data . . . .')
fid = fopen(fname, 'rb');
data = fread(fid, [xdim * ydim * zdim * timepoints], type);
disp(['max: ' num2str(max(data(:))) ', min: '  num2str(min(data(:))) ])
fclose(fid);
disp('BRIK2MGH: done reading raw data !')
disp(' ')

disp('BRIK2MGH: reshaping raw data to 3D+time matrix . . . .')
V = reshape(data,xdim,ydim,zdim,timepoints);
disp('BRIK2MGH: done reshaping raw data to 3D+time matrix !')
disp(' ')

disp('BRIK2MGH: writing MGH-style image files . . . .')
disp(' ')

bfile =  zeros(xdim,ydim,timepoints);

% cycle over slices and make a volume
% where the new zdimension is timpoints

for i = 1:zdim

 % make file of timepoints
 bfile(:,:,:) = V(:,:,i,:);

 % create the bfile name
 suffix = ['.b' type];
 outfname = sprintf('%s%03d%s%s',foutstem,i-1,suffix);

 % write out the bfile
 
 disp(['BRIK2MGH: writing file ' num2str(outfname) ', max: ' num2str(max(bfile(:))) ', min: ' num2str(min(bfile(:))) ])
 fid = fopen(outfname, 'w');
 fwrite(fid, bfile(:), type);
 fclose(fid);

 % write out a header file
 hdrfname = sprintf('%s%03d%s',foutstem,i-1,'.hdr');
 fid = fopen(hdrfname,'w');
 fprintf(fid,'%d %d %d 0\n',xdim, ydim, timepoints);
 fclose(fid);

end

disp(' ')
disp('Done!')


      

	
      
        
         
