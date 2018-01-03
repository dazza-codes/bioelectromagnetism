function fast_fwhm(input_files,sd_file,fwhm_file)

n=size(input_files,1);
d=fmris_read_image(input_files(1,:),0,0);
m=0;
for i=1:n
   d=fmris_read_image(input_files(i,:),1:d.dim(3),1);
   m=m+d.data;
end
m=m/n;

v0=0;
for i=1:n
   d=fmris_read_image(input_files(i,:),1:d.dim(3),1);
   v0=v0+(d.data-m).^2;
end
out.data=sqrt(v0/(n-1))*sqrt(2/100);
out.file_name=sd_file;
out.dim=[d.dim(1:3) 0];
out.parent_file=input_files(1,:);
fmris_write_image(out,1:d.dim(3),1);

v00=1./sqrt(v0+(v0<=0)).*(v0>0);
v11=0;
v12=0;
v13=0;
v22=0;
v23=0;
v33=0;
for i=1:n
   i
   d=fmris_read_image(input_files(i,:),1:d.dim(3),1);
   [d1,d2,d3]=gradient((d.data-m).*v00);
   v11=v11+d1.*d1;
   v12=v12+d1.*d2;
   v13=v13+d1.*d3;
   v22=v22+d2.*d2;
   v23=v23+d2.*d3;
   v33=v33+d3.*d3;
end
detlam=v11.*(v22.*v33-v23.^2)-v12.*(v12.*v33-v23.*v13)+v13.*(v12.*v23-v22.*v13);
out.data=sqrt(4*log(2)./(detlam+(detlam<=0)).^(1/3))*prod(abs(d.vox))^(1/3).*(detlam>0); 
out.file_name=fwhm_file;
out.dim=[d.dim(1:3) 2];
fmris_write_image(out,1:d.dim(3),1);
out.data=sqrt(detlam+(detlam<=0)).*(detlam>0);
fmris_write_image(out,1:d.dim(3),2);

return
