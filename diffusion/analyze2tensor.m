fname1='/local_mount/space/tensor/1/users/mette/MS/HUTTENBAUER3403918/3403918_10898_002_image1.ah/temp2.bshort';
fname2='/local_mount/space/tensor/1/users/mette/MS/HUTTENBAUER3403918/3403918_10898_002_image1.ah/temp2_tensor.bfloat';

b=1000;


fid=fopen(fname1,'r');
buff=fread(fid,inf,'uint16');
fclose(fid);

nslc=18;
a=squeeze(reshape(buff,[128 128 7 nslc]));

for k=11:11 %nslc
    sprintf('calculating diffusion coefficients for slice #%d',k)
    for i=1:128
	for j=1:128;
	    for l=1:7
		if(isinf(a(i,j,l,k))|isnan(a(i,j,l,k))) a(i,j,l,k)=0; 
		end
	    end
            for l=1:6 
              d(i,j,l,k)=-1/b*log(a(i,j,l,k)/a(i,j,7,k));
            end
            dd(i,j,1,k)=-1/(4*b)*(log(a(i,j,1,k))+log(a(i,j,2,k))-log(a(i,j,3,k))-log(a(i,j,4,k))+log(a(i,j,5,k))+log(a(i,j,6,k))-2*log(a(i,j,7,k)));
            dd(i,j,2,k)=-1/(4*b)*(log(a(i,j,1,k))-log(a(i,j,2,k)));
            dd(i,j,3,k)=-1/(4*b)*(-log(a(i,j,5,k))+log(a(i,j,6,k)));
            dd(i,j,4,k)=-1/(4*b)*(log(a(i,j,1,k))+log(a(i,j,2,k))+log(a(i,j,3,k))+log(a(i,j,4,k))-log(a(i,j,5,k))-log(a(i,j,6,k))-2*log(a(i,j,7,k)));
            dd(i,j,5,k)=-1/(4*b)*(-log(a(i,j,3,k))+log(a(i,j,4,k)));
            dd(i,j,6,k)=-1/(4*b)*(-log(a(i,j,1,k))-log(a(i,j,2,k))+log(a(i,j,3,k))+log(a(i,j,4,k))+log(a(i,j,5,k))+log(a(i,j,6,k))-2*log(a(i,j,7,k)));

	    tmpx(i,j,:,k)=[dd(i,j,1,k) dd(i,j,2,k) dd(i,j,3,k) dd(i,j,4,k) dd(i,j,5,k) dd(i,j,6,k)];
	    
	end
    end
end
tmpx=reshape(tmpx,[128 128 6*11]);

fid=fopen(fname2,'w');
buff=fwrite(fid,tmpx,'float32');
fclose(fid);

clear all




