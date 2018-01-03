function rgbdata=spm_diffusion_suconvert
% 
hdr_size=3952;  %size of GE hdr in int16 units
  
t2_ims = spm_get([1,100],'TENSOR_B0*.0*','select low-b images');
outp_files =  spm_get([1,100],'*.outp*','select .outp files');

% read in images

[path t2_stem foo foo2]=fileparts(t2_ims(1,:));

fa_ims = spm_get('files',path,'TENSOR_FA*.0*');
dw_ims = spm_get('files',path,'TENSOR_DW*.0*');
ci_ims = spm_get('files',path,'tensor_CI*.0*');
tr_ims = spm_get('files',path,'TENSOR_TR*.0*');

[path fa_stem foo foo2]=fileparts(fa_ims(1,:));
[path dw_stem foo foo2]=fileparts(dw_ims(1,:));
[path ci_stem foo foo2]=fileparts(ci_ims(1,:));
[path tr_stem foo foo2]=fileparts(tr_ims(1,:));

[su_hdr,ex_hdr,se_hdr,im_hdr,pix_hdr,im_offset] = GE_readHeader(t2_ims(1,:));
img_size=[im_hdr.imatrix_X,im_hdr.imatrix_Y];
hdr_size=pix_hdr.img_hdr_length/2;
voxdims=[im_hdr.pixsize_X,im_hdr.pixsize_Y,im_hdr.slthick];

% data structure containing all the data
%
% (1) - t2 image
% (2) - FA image
% (3) - DWI image
% (4) - lattice anisotropy image
% (5) - trace image
% (6:8) - RGB image for 1st eigenvector
% (9:11) - 
data=zeros(img_size(1),img_size(2),size(t2_ims,1),5);  % 128 x 128 x nslices x 2 (t2,fa)
rgbdata=zeros(img_size(1),img_size(2),size(t2_ims,1),3);
xydata=zeros(img_size(1),img_size(2),size(t2_ims,1),3);
fadata=zeros(img_size(1),img_size(2),size(t2_ims,1));

% read in data from .outp files

[outp_path outp_stem foo foo2]=fileparts(outp_files(1,:));

for f=1:size(t2_ims,1),
  % first, read in t2 image
   fid=fopen(t2_ims(f,:),'r','b');
   imghdr=fread(fid,hdr_size,'int16');    
   data(:,:,f,1)=fread(fid,img_size,'int16'); 
   fclose(fid);
  % then, read in DWI image 
   fid=fopen(dw_ims(f,:),'r','b');
   imghdr=fread(fid,hdr_size,'int16');    
   data(:,:,f,3)=fread(fid,img_size,'int16'); 
   fclose(fid);
  % then, read in CI image 
   fid=fopen(ci_ims(f,:),'r','b');
   imghdr=fread(fid,hdr_size,'int16');    
   data(:,:,f,4)=fread(fid,img_size,'int16'); 
   fclose(fid);
  % then, read in TR image 
   fid=fopen(tr_ims(f,:),'r','b');
   imghdr=fread(fid,hdr_size,'int16');    
   data(:,:,f,5)=fread(fid,img_size,'int16');  
   fclose(fid);

  % get RGB data from outp file

   [outp_path outp_stem foo foo2]=fileparts(outp_files(f,:));
   [Dot X Y Z Emax Emin FA RGB] = spm_diffusion_loaddot(outp_stem,0,img_size(1));
   rgbdata(:,:,f,1)=squeeze(RGB(:,:,1));
   rgbdata(:,:,f,2)=squeeze(RGB(:,:,2));
   rgbdata(:,:,f,3)=squeeze(RGB(:,:,3));
   xydata(:,:,f,1)=X;
   xydata(:,:,f,2)=Y;
   xydata(:,:,f,3)=Z;
   data(:,:,f,2)=FA;
   
end;

vox=voxdims;
dim=[img_size(1),img_size(2),size(t2_ims,1),4];

[outp_path outp_stem foo foo2]=fileparts(outp_files(1,:));

V=struct(...
         'fname',  sprintf('%s.img',t2_stem),...
	 'dim',  dim,...
	 'mat', [-1*vox(1) 0 0 (vox(1)*dim(1))*0.5; 0 -1*vox(2) 0 (vox(2)*dim(2))*0.45; 0 0 vox(3) -1*(vox(3)*dim(3))*0.33; 0 0 0 1],...
	 'pinfo',[1;0;0]...
);

spm_write_vol(V,squeeze(data(:,:,:,1)));

V.fname=sprintf('%s.img',dw_stem);
spm_write_vol(V,squeeze(data(:,:,:,3)));

V.fname=sprintf('%s.img',ci_stem);
spm_write_vol(V,squeeze(data(:,:,:,4)));

V.fname=sprintf('%s.img',tr_stem);
spm_write_vol(V,squeeze(data(:,:,:,5)));

V.dim=[img_size(1),img_size(2),size(t2_ims,1),64];

V.fname=sprintf('%s.img',fa_stem);
spm_write_vol(V,squeeze(data(:,:,:,2)));

V.fname=sprintf('%s_col1_1.img',outp_stem);
spm_write_vol(V,squeeze(rgbdata(:,:,:,1)));

V.fname=sprintf('%s_col1_2.img',outp_stem);
spm_write_vol(V,squeeze(rgbdata(:,:,:,2)));

V.fname=sprintf('%s_col1_3.img',outp_stem);
spm_write_vol(V,squeeze(rgbdata(:,:,:,3)));


V.fname=sprintf('%s_X.img',outp_stem);
spm_write_vol(V,squeeze(xydata(:,:,:,1)));

V.fname=sprintf('%s_Y.img',outp_stem);
spm_write_vol(V,squeeze(xydata(:,:,:,2)));

V.fname=sprintf('%s_Z.img',outp_stem);
spm_write_vol(V,squeeze(xydata(:,:,:,3)));


