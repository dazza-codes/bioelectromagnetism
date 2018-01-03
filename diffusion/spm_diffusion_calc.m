function spm_diffusion_calc(images,b)
% spm_diffusion_calc: part of the SPM diffusion tensor toolbox
% Russ Poldrack, 11/19/00
%
% - This routine calculates the diffusion tensor and fractional 
%   anisotropy based upon a set of 6 diffusion-weighted images along 
%   with a non-diffusion-weighted (low-b) image.  The calculation is based
%   on code provided by Mette Weigell at the MGH-NMR Center.
%
% The particular gradient directions are hard-coded into the calculation
% (as listed below).  If someone would like to create a general version
% of this calculation code such that gradient directions could be specified
% on the fly, the code would become much more generally useful. 
%
% The calculations do not include any corrections for gradient cross-terms,
% and no corrections are applied for image distortion resulting from 
% eddy currents.  For a properly designed (balanced refocused) pulse sequence 
% this shouldn't be a problem.
%
% Arguments:
% images: this should specify a set of 7 images, with the following 
% diffusion weightings (images must in this specific order):
%
%     X   Y   Z
% 1 - 0   0   0 (i.e., no diffusion weighting)
% 2 - 1   1   0
% 3 - 1  -1   0
% 4 - 0   1   1
% 5 - 0  -1   1
% 6 - 1   0   1
% 7 - -1  0   1
%   
% b: b value (i.e., amount of diffusion weighting)
%
% Effects:
% creates the following files:
%  fa.img - fractional anisotropy
%  col1_[1-3].img - images representing the R,G,and B components respectively
%   of a diffusion colormap for the first eigenvector of the tensor
%  col3_[1-3].img - images representing the R,G,and B components respectively
%   of a diffusion colormap for the third eigenvector of the tensor
%  eig1_{X,Y,Z}.img - images representing the X,Y, and Z components of the 
%   first eigenvector respectively.
%  eig.bfloat - binary image containing tensor and eigenvector information
%   (for use with in-house MGH tools)
%
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% The SPM diffusion toolbox is distributed freely under the GNU Public License.
%

vols=spm_vol(images);
[data XYZ]=spm_read_vols(vols);
a=zeros(size(data));
a(:,:,:,1:6)=data(:,:,:,2:7);
a(:,:,:,7)=data(:,:,:,1);  % put low-b data into 7

% set up b matrix
if isempty('b'),b=600;end;


% $$$ % pick an example voxel: 71,49,12  (genu of CC)
% $$$      S=zeros(6,1);
% $$$      location=[71 49 12];
% $$$      S_lowb=data(location(1),location(2),location(3),1);
% $$$      for wtg=1:6,
% $$$ 	S(wtg)=data(location(1),location(2),location(3),wtg+1);
% $$$      end;
% $$$      fprintf('%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n',...
% $$$ 	     S_lowb,S);


% a=squeeze(reshape(buff,[128 128 7 nslc]));
[xsize,ysize,nslc,foo]=size(a);

% note - need to switch k and l indices

%  % FOR DEBUGGING ONLY
tmpx=zeros(xsize,ysize,nslc,6);
lowb_thresh=50;

spm_progress_bar('Init',nslc*2);
spm('Pointer','Watch');
prog_counter=0;

for k=1:nslc
    fprintf('calculating diffusion coefficients for slice #%d\n',k);
%    imagesc(a(:,:,k,1));
    for i=1:xsize
	for j=1:ysize;
          warning off;
	  for l=1:7
		if(isinf(a(i,j,k,l))|isnan(a(i,j,k,l))) a(i,j,k,l)=0; end
	  end;
            if a(i,j,k,7) > lowb_thresh
              for l=1:6 
                d(i,j,k,l)=-1/b*log(a(i,j,k,l)/a(i,j,k,7));
              end
              dd(i,j,k,1)=-1/(4*b)*(log(a(i,j,k,1))+log(a(i,j,k,2))-log(a(i,j,k,3))-log(a(i,j,k,4))+log(a(i,j,k,5))+log(a(i,j,k,6))-2*log(a(i,j,k,7)));
              dd(i,j,k,2)=-1/(4*b)*(log(a(i,j,k,1))-log(a(i,j,k,2)));
              dd(i,j,k,3)=-1/(4*b)*(-log(a(i,j,k,5))+log(a(i,j,k,6)));
              dd(i,j,k,4)=-1/(4*b)*(log(a(i,j,k,1))+log(a(i,j,k,2))+log(a(i,j,k,3))+log(a(i,j,k,4))-log(a(i,j,k,5))-log(a(i,j,k,6))-2*log(a(i,j,k,7)));
              dd(i,j,k,5)=-1/(4*b)*(-log(a(i,j,k,3))+log(a(i,j,k,4)));
              dd(i,j,k,6)=-1/(4*b)*(-log(a(i,j,k,1))-log(a(i,j,k,2))+log(a(i,j,k,3))+log(a(i,j,k,4))+log(a(i,j,k,5))+log(a(i,j,k,6))-2*log(a(i,j,k,7)));

	      tmpx(i,j,k,:)=[dd(i,j,k,1) dd(i,j,k,2) dd(i,j,k,3) dd(i,j,k,4) dd(i,j,k,5) dd(i,j,k,6)];
	    else
              for l=1:6 
                d(i,j,k,l)=0;
              end
	        dd(i,j,k,1)=0;
	        dd(i,j,k,2)=0;
	        dd(i,j,k,3)=0;
	        dd(i,j,k,4)=0;
	        dd(i,j,k,5)=0;
	        dd(i,j,k,6)=0;

	      tmpx(i,j,k,:)=[dd(i,j,k,1) dd(i,j,k,2) dd(i,j,k,3) dd(i,j,k,4) dd(i,j,k,5) dd(i,j,k,6)];     
	    end;
	   warning on;
	end
    end
    prog_counter=prog_counter+1;
    spm_progress_bar('Set',prog_counter);
end

% code adapted from tensor2eig.m

clear a;
a=tmpx;
clear tmpx tmpy;

for k=1:nslc
    fprintf('calculating eigensystem and FA for slice #%d\n',k);
    for i=1:xsize
	for j=1:ysize
	    for l=1:6
		if(isinf(a(i,j,k,l))|isnan(a(i,j,k,l))) a(i,j,k,l)=0; 
		end
	    end

	    if isempty(find(a(i,j,k,:))), 
	         D=[1 0 0; 0 1 0; 0 0 1]; 
            else
	      D=[a(i,j,k,1) a(i,j,k,2) a(i,j,k,3);
	       a(i,j,k,2) a(i,j,k,4) a(i,j,k,5);
	       a(i,j,k,3) a(i,j,k,5) a(i,j,k,6)];
            end;
	    
	    [v,d]=eig(D);
	    d=diag(d);
	    dq=flipud(sort(d));
	    for m=1:3
	        for n=1:3 
		    if(d(n)==dq(m))
		        qz(m)=n;
		    end
		end
	    end

	    d=d(qz);    
	    vtemp=zeros(3,3);
	    for m=1:3
		vtemp(:,m)=v(:,qz(m));
	    end
	    v=vtemp;
	    tmpx(i,j,k,:)=[d' v(:,1)' v(:,2)' v(:,3)'];
	    if(sum(d)==0) 
		FA=0;
	    else
		L=d;
		FA=sqrt(3/2)*sqrt(1/3*((L(1)-L(2))^2+(L(2)-L(3))^2+(L(3)-L(1))^2))/sqrt(L(1)^2+L(2)^2+L(3)^2);
	    end
	tmpy(i,j,k,:)=FA;
	end
    end
    prog_counter=prog_counter+1;
    spm_progress_bar('Set',prog_counter);
end

tmpy=reshape(tmpy,[xsize ysize nslc]);

tmpx=reshape(tmpx,[xsize ysize 12*nslc]);

DIM=vols(1).dim;
VOX=[vols(1).mat(1,1)  vols(1).mat(2,2)  vols(1).mat(3,3)];
ORIGIN = (vols(1).mat\[0 0 0 1]')'

spm_hwrite('fa.img',DIM,VOX,1,16,0,ORIGIN,'tensor');
fid=fopen('fa.img','w'); fwrite(fid,tmpy,'float'); fclose(fid);


fid=fopen('eig.bfloat','w');
buff=fwrite(fid,tmpx,'float32');
fclose(fid);

% now write color maps - each in three files

buff1=tmpx;
buff2=tmpy;
buff3=squeeze(data(:,:,:,1));  % low-b
buff4=mean(data(:,:,:,2:7),4);



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

spm_hwrite('col3_1.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');
spm_hwrite('col3_2.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');
spm_hwrite('col3_3.img',DIM,VOX,1,16,0,ORIGIN,'tensor colormap');

fid=fopen('eig1_X.img','w'); fwrite(fid,E1x,'float'); fclose(fid);
fid=fopen('eig1_Y.img','w'); fwrite(fid,E1y,'float'); fclose(fid);
fid=fopen('eig1_Z.img','w'); fwrite(fid,E1z,'float'); fclose(fid);

spm_hwrite('eig1_X.img',DIM,VOX,1,16,0,ORIGIN,'tensor eig1 x');
spm_hwrite('eig1_Y.img',DIM,VOX,1,16,0,ORIGIN,'tensor eig1 y');
spm_hwrite('eig1_Z.img',DIM,VOX,1,16,0,ORIGIN,'tensor eig1 z');

tmpvol1=zeros(DIM(1:3));
tmpvol2=zeros(DIM(1:3));
tmpvol3=zeros(DIM(1:3));

col3=cat(4,abs(E3x).*sqrt(Df2scaled).*mask,abs(E3y).*sqrt(Df2scaled).*mask,abs(E3z).*sqrt(Df2scaled).*mask);

for x=1:DIM(1),
  for y=1:DIM(2),
    for z=1:DIM(3),
      tmpvol1(x,y,z)=col3(y,x,z,1);
      tmpvol2(x,y,z)=col3(y,x,z,2);
      tmpvol3(x,y,z)=col3(y,x,z,3);
    end;
  end;
end;
fid=fopen('col3_1.img','w'); fwrite(fid,tmpvol1,'float'); fclose(fid);
fid=fopen('col3_2.img','w'); fwrite(fid,tmpvol2,'float'); fclose(fid);
fid=fopen('col3_3.img','w'); fwrite(fid,tmpvol3,'float'); fclose(fid);





spm_progress_bar('Clear');
spm('Pointer','Arrow');



