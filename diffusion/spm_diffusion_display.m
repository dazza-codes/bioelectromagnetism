	
%img   = spm_get(1,'fa.img','select fa.img');
img='fa.img';	
v=spm_vol(img);
[data xyz]=spm_read_vols(v);
whos
x=[8];
xyz=[0;0;0];
clf;
%hMIPax = spm_mip_ui(x,xyz,v.mat,v.dim(1:3)',gcf);

H = spm_orthviews('Image','fa.img',[0.0 0.45 1 0.55])

%H = spm_diffusion_orthviews('Image','fa.img',[0.0 0.0 1 0.55])






















