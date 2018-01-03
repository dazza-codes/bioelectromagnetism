function [FA,FAscaled,trace,col1,col3,lowb,DWI]=load_tensor_data(eigfile,avgdwi)
	global FA FAscaled trace col1 col3 lowb DWI;

	[dname eigfile ext]=dirname(eigfile);
	fprintf('Loading diffusion data...\n');

fid=fopen(sprintf('%s/eig.bfloat',dname),'r');
buff1=fread(fid,inf,'float32');
fclose(fid);

vols=spm_vol(sprintf('%s/fa.img',dname));
[buff2 XYZ]=spm_read_vols(vols);
[xsize,ysize,nslc]=size(buff2);
fa_vol=vols;

vols=spm_vol(sprintf('%s/%s',dname,avgdwi));
[buff3 XYZ]=spm_read_vols(vols);

vols=spm_vol(sprintf('%s/%s',dname,avgdwi));
[buff4 XYZ]=spm_read_vols(vols);


imdeck=squeeze(reshape(buff1,[xsize ysize 12*nslc]));
imdeck=permute(imdeck,[2,1,3]);
[npix,nim,imd]=getimd(imdeck);

imdeck2=squeeze(reshape(buff2,[xsize ysize nslc]));
imdeck2=permute(imdeck2,[2,1,3]);
[npix,nim,imd]=getimd(imdeck2);

imdeck3=squeeze(reshape(buff3,[xsize ysize  nslc]));
imdeck3=permute(imdeck3,[2,1,3]);
[npix,nim,imd]=getimd(imdeck3);

imdeck4=squeeze(reshape(buff4,[xsize ysize nslc]));
imdeck4=permute(imdeck4,[2,1,3]);
[npix,nim,imd]=getimd(imdeck3);

npar=12;
n=nslc*npar;
npar2=1;
n2=nslc*npar2;

L1=imdeck(:,:,1:nslc);
L2=imdeck(:,:,nslc+1:2*nslc);
L3=imdeck(:,:,2*nslc+1:3*nslc);

E1x=imdeck(:,:,3*nslc+1:4*nslc);
E1y=imdeck(:,:,4*nslc+1:5*nslc);
E1z=imdeck(:,:,5*nslc+1:6*nslc);

E3x=imdeck(:,:,6*nslc+1:7*nslc);
E3y=imdeck(:,:,7*nslc+1:8*nslc);
E3z=imdeck(:,:,8*nslc+1:9*nslc);

FA=imdeck2(:,:,:);
lowb=imdeck3(:,:,:);
DWI=imdeck4(:,:,:);

Df2=abs(L2-L3);

FAmax=max(reshape(FA,npix,nslc));
FAmax=reshape(repmat(FAmax,npix,1),[imd(1) imd(2) nslc]);
FAscaled=FA./FAmax;
trace=abs(L1+L2+L3);

mask=((lowb>10) & (FA<1) & (trace<3.5*10^-3));	

Df2max=max(reshape(Df2,npix,nslc));
Df2max=reshape(repmat(Df2max,npix,1),[imd(1) imd(2) nslc]);
Df2scaled=sqrt(Df2./Df2max);

eigval=cat(4,abs(L1),abs(L2),abs(L3));

col1=cat(4,abs(E1x).*sqrt(FAscaled).*mask,abs(E1y).*sqrt(FAscaled).*mask,abs(E1z).*sqrt(FAscaled).*mask);

DIM=fa_vol(1).dim;
VOX=[fa_vol(1).mat(1,1)  fa_vol(1).mat(2,2)  fa_vol(1).mat(3,3)];
ORIGIN=[DIM(1)/2  DIM(2)/2  DIM(3)/2]

spm_hwrite('col1_1.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');
spm_hwrite('col1_2.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');
spm_hwrite('col1_3.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');

tmpvol1=zeros(DIM(1:3));
tmpvol2=zeros(DIM(1:3));
tmpvol3=zeros(DIM(1:3));

for x=1:DIM(1),
  for y=1:DIM(2),
    for z=1:DIM(3),
      tmpvol1(x,y,z)=col1(y,x,z,1);
      tmpvol2(x,y,z)=col1(y,x,z,2);
      tmpvol3(x,y,z)=col1(y,x,z,3);
    end;
  end;
end;
fid=fopen('col1_1.img','w'); fwrite(fid,tmpvol1,'float'); fclose(fid);
fid=fopen('col1_2.img','w'); fwrite(fid,tmpvol2,'float'); fclose(fid);
fid=fopen('col1_3.img','w'); fwrite(fid,tmpvol3,'float'); fclose(fid);


col3=cat(4,abs(E3x).*sqrt(Df2scaled).*mask,abs(E3y).*sqrt(Df2scaled).*mask,abs(E3z).*sqrt(Df2scaled).*mask);


clear buff* imdeck* *sizeL* E*  Df* FAmax imd tmp* n* imdir* fn* ans fid mask
whos
rip

















