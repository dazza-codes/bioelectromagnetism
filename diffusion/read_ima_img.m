function SpmImgStruct=read_ima_img(file)
% 
% Russ Poldrack, 8/15/2000
% 



     ImaHdr=read_siemens_header(file);

        DIM(1) = ImaHdr.h_G19_Acq3_Mr_BaseRawMatrixSize; %matrix size of each image 
        DIM(2) = DIM(1); 
        DIM(3) = 1; % num slices unknown 
        slvox=ImaHdr.h_G21_Rel1_CM_FoV_Width/ImaHdr.h_G19_Acq3_Mr_BaseRawMatrixSize;

        VOX(1) = ImaHdr.h_G21_Rel1_CM_FoV_Width/ImaHdr.h_G19_Acq3_Mr_BaseRawMatrixSize;
        VOX(2) = VOX(1); 
        VOX(3) = ImaHdr.h_G18_Acq_SliceThickness; 
         
        % ORIGIN is just the center of the image matrix 
        % more fancy ways of doing this exists 

        z = NaN; % cannot calculate because slice thickness is unknown 
        x = ImaHdr.h_G19_Acq3_Mr_BaseRawMatrixSize*0.5; 
        y = ImaHdr.h_G19_Acq3_Mr_BaseRawMatrixSize*0.5; 

        ORIGIN = [x ; y; z]; 

        % read the image 
        fid = fopen(file,'r','b'); % 'b' argument denotes big-endian (from siemens/sun)
        fseek(fid,6144,'bof'); 
        V = fread(fid, [DIM(1) DIM(2)], 'int16'); 

        fclose(fid); 




SpmImgStruct=struct( ...
       'ImaHdr',     ImaHdr, ...
       'ORIGIN',     ORIGIN, ...
       'VOX',        VOX, ...
       'DIM',        DIM, ...
       'V',          V);







