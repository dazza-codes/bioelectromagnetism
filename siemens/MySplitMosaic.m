 function SplitMosaic(inpath, outpath, uid, start, stop, nslice) 
% 
% splits Siemens Mosaic images and writes Analyze images 
% if outpath does not exist it is created 
% 
% 990525 GH 
% Modified by Eman Ghobrial Feb 2000
% 
% calls the following subroutines: 
% ZeroPad(prefix,npos,fileno) 
% spm_hwrite(P,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP) 
% 

 if nargin < 7 
    inpath = input('Siemens Mosaic path: ','s'); 
    outpath= input('Output path : ','s'); 
    if ~exist(outpath) 
      unix(['mkdir ' outpath]); 
    end; 
    uid = input('Patient no: '); 
    start = input('First series no: '); 
    startno = num2str(start); 

    stop = input('Last series no: '); 
    nslice = input('number of slices: '); 
  end 
  
  
  voltot = 0; 
  for series = start:stop 
     voltot = voltot + 1;
     
         % First image no: 
         nchar = length(sprintf('%d-%d-',uid,series))+1;
         imageList = [];
       if strcmp(computer,'PCWIN')
          d = dir(sprintf('%s\\%s-%d-*.ima', inpath, num2str(uid), series));
          numofima = length(d)
       for i=1:length(d)
          filename = d(i).name;
			 imanum(i)=str2num(filename(nchar:length(filename)-4));
		 end
       %first image number
       imageList = sort(imanum);
       firstima = imageList(1)
%	firstima = str2num(d(1).name(nchar:length(d(1).name)-4))
    else
     		d = dir(sprintf('%s/%s-%d-*.ima', inpath, num2str(uid), series)) 
       numofima = length(d);
       for i=1:length(d)
          filename = d(i).name;
			 imanum(i)=str2num(filename(nchar:length(filename)-4));
		 end
	%first image number
		imageList = sort(imanum);
      firstima = imageList(1)

    end    

%go through the images     
  	for count = 1:numofima   
    if strcmp(computer,'PCWIN')
   %file = [inpath '\\' num2str(uid) '-' num2str(series) '-' imano '.ima']; 
	 file = sprintf('%s\\%s-%s-%d.ima',inpath,num2str(uid),num2str(series),imageList(count));
 	else
    	 file = sprintf('%s/%s-%s-%d.ima',inpath,num2str(uid),num2str(series),imageList(count));
	 end    

   if strcmp(computer,'PCWIN')      
      fid = fopen(file,'r','b');
   else   
      fid = fopen(file,'r'); 
   end   
    if fid < 0 
      fprintf('File %s not found\n',file); 
		count = count+ 1; 
    else 
        fprintf('Opening file %s \n',file); 
        % Get header information, 
        fseek(fid,2864,'bof'); 
        mtrx = fread(fid,1,'int32'); 
        fseek(fid,3744,'bof'); 
        fov = fread(fid,1,'double'); 

        slvox = fov/mtrx; 

        fseek(fid,5000,'bof'); 
        imgvox = fread(fid,1,'double'); 
        ncol = slvox/imgvox; 
        nrow = ncol; 

        DIM(1) = mtrx; %matrix size of each image 
        DIM(2) = DIM(1); 
        DIM(3) = nslice; 

        VOX(1) = slvox; 
        VOX(2) = VOX(1); 
        fseek(fid,1544,'bof'); 
        VOX(3) = fread(fid,1,'double'); 

         
        % ORIGIN is just the center of the image matrix 
        % more fancy ways of doing this exists 

        z = nslice/2; 
        x = mtrx*0.5; 
        y = mtrx*0.5; 

        ORIGIN = [x ; y; z]; 

        % read the image 
        fseek(fid,6144,'bof'); 
        V = fread(fid, [ncol*mtrx nrow*mtrx], 'int16'); 

        % mesh(V); 
        % view(2); 

        for j = 1:nrow 
          for k = 1:ncol 
            slice = (j-1)*ncol + k; 
            if slice <= nslice 
              A=V(k*mtrx-mtrx+1:k*mtrx, j*mtrx-mtrx+1:j*mtrx); 
              img(:,:,nslice-slice+1)=flipud(fliplr(A)); 
              %figure; 
              %string = ['this is sliceno: ' num2str(slice)]; 
              %mesh(img(:,:,slice)), title(string); 
              %view(2); 
            end; 
          end; 
        end; 
        fclose(fid); 
        if strcmp(computer,'PCWIN')
            prefix =[outpath '\\' num2str(uid) ZeroPad('-',4,series) ZeroPad('-',4,imageList(count))]; 
			else 
            prefix =[outpath '/' num2str(uid) ZeroPad('-',4,series) ZeroPad('-',4,imageList(count))]; 
         end   
        spm_hwrite(prefix, DIM, VOX, 8, 4, 0, ORIGIN, 'Siemens Vision'); 
        fid = fopen([prefix '.img'],'w'); 
        if (fid > 0) 
                % fprintf('Writing Analyze file no %s \n',prefix); 
                fwrite(fid, img, 'int16'); 
                fclose(fid); 
        end 
     end; %if file was opened 
     
             count = count +1; 

   end; %for each ima  
  end; %for each series 

vol = voltot; 
vox = VOX; 
dim = DIM(1); 
fprintf('\n'); 
fprintf('*********************************************************\n'); 
fprintf('Total no of volumes is: %5.0f\n',vol); 
fprintf('\n'); 
fprintf('Voxelsize, x direction: %2.2f, y: %2.2f, z:%2.2f\n',vox); 
fprintf('\n'); 
fprintf('Matrixsize: %3.0f\n',dim); 
fprintf('*********************************************************\n'); 
fprintf('\n'); 

return
