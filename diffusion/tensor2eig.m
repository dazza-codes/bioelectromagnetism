%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  filename = fname1 for diffusion tensors
%  filename = fname2 for eigenvalues and eigenvectors
%  filename = fname3 for FA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imdir_inp='/jaz1/users/mette/GBM/3737163_2/data.1895.5.1/run';
oldornew='4';
%% 2 for old data, 4 for new data or S for Siemens data

oldornew=sprintf('%s',oldornew);
fname1=sprintf('%s/temp%s_tensor.bfloat',imdir_inp,oldornew);
fname2=sprintf('%s/temp%seig_1.bfloat',imdir_inp,oldornew);
fname3=sprintf('%s/temp%sfa_1.bfloat',imdir_inp,oldornew);

fid=fopen(fname1,'r');
buff=fread(fid,inf,'float32');
fclose(fid);

if (oldornew=='2')
    %nslc=22;
    nslc=18;
    else
    if (oldornew=='4')
	nslc=23;
	%nslc=21;
	b=1000;
	else
	nslc=28;
	b=700;
    end
end;

a=squeeze(reshape(buff,[128 128 6 nslc]));
for k=1:nslc
    sprintf('calculating eigensystem and FA for slice #%d',k)
    for i=1:128
	for j=1:128;
	    for l=1:6
		if(isinf(a(i,j,l,k))|isnan(a(i,j,l,k))) a(i,j,l,k)=0; 
		end
	    end

	    D=[a(i,j,1,k) a(i,j,2,k) a(i,j,3,k);
	       a(i,j,2,k) a(i,j,4,k) a(i,j,5,k);
	       a(i,j,3,k) a(i,j,5,k) a(i,j,6,k)];

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
	    tmpx(i,j,:,k)=[d' v(:,1)' v(:,2)' v(:,3)'];
	    if(sum(d)==0) 
		FA=0;
	    else
		L=d;
		FA=sqrt(3/2)*sqrt(1/3*((L(1)-L(2))^2+(L(2)-L(3))^2+(L(3)-L(1))^2))/sqrt(L(1)^2+L(2)^2+L(3)^2);
	    end
	tmpy(i,j,:,k)=FA;
	end
    end
end
tmpx=reshape(tmpx,[128 128 12*nslc]);
tmpy=reshape(tmpy,[128 128 nslc]);

fid=fopen(,'w');
buff=fwrite(fid,tmpx,'float32');
fclose(fid);

fid=fopen(fname3,'w');
buff=fwrite(fid,tmpy,'float32');
fclose(fid);

clear all



