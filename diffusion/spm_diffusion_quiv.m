function spm_diffusion_quiv(anat_img,fa_img,x_img,y_img,z_img)
% spm_diffusion_quiv: part of the SPM diffusion toolbox
% Russ Poldrack, 11/19/00
%
% - This routine creates a vector map of diffusion eigenvectors
%   overlaid on an anatomical image.  In its current state, it is
%   mildly interactive, in the sense that different slices can be 
%   selected by specifying the slice number.  In the future I hope
%   to integrate it with the SPM functions that allow interactive
%   specification of location, and perhaps also with the color map
%   display.
%
% Arguments:
% anat_img: anatomical image
% fa_img: FA image
% {x,y,z}_img: eigenvector images for X, Y, and Z directions
%
% Effects:
% Displays vector map overlaid on anatomical image
%
% CVS repository and development information can be found at:
% http://spm-toolbox.sourceforge.net 
%
% The SPM diffusion toolbox is distributed freely under the GNU Public License.
%

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Display',0);
SPMid = spm('FnBanner',mfilename,'2.15');
spm_help('!ContextHelp',[mfilename,'.m']);

%anat_img=spm_get(1,'.img','select an anatomical image');
%fa_img=spm_get(1,'.img','select an FA image');
%xyz_img=spm_get(3,'.img','select [X,Y,Z] images');

xyz_img=[x_img;y_img;z_img];
anat_v=spm_vol(anat_img);
xyz_v=spm_vol(xyz_img);
fa_v=spm_vol(fa_img);
anat=spm_read_vols(anat_v);
xyz=spm_read_vols(xyz_v);
xya=abs(xyz);
fa=spm_read_vols(fa_v);
fa_min=0.25;


slice=1;
slice=spm_input('slice number to view (0 to quit)',anat_v.dim(3));
while slice > 0,
  figure(Fgraph);
  spm_figure('Clear',Fgraph);
  hold on;
  axis square;
  axis off;
  rotimg=rot90(anat(:,:,slice),-1);
rotimg=fliplr(rotimg);
  imagesc(rotimg);
  quiv_data=zeros(4,anat_v.dim(1)*anat_v.dim(2));

 % threshold using FA
  counter=1;
  for x=1:anat_v.dim(1),
    for y=1:anat_v.dim(2),
      if fa(x,y,slice)>fa_min, 
          quiv_data(:,counter)=[x;y;xyz(x,y,slice,1)*fa(x,y,slice);xyz(x,y,slice,2)*fa(x,y,slice)];
          counter=counter+1;
      end;
    end;
  end;
quiv_data=quiv_data(:,1:counter);
  quiver(quiv_data(1,:),quiv_data(2,:),quiv_data(3,:),quiv_data(4,:),'r.-');
  hold off;
  slice=spm_input('slice number to view (0 to quit)',anat_v.dim(3));

end;

